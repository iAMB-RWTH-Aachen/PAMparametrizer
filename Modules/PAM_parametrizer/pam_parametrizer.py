from typing import Union, Tuple, Iterable
import numpy as np
import os
import time
import random
import matplotlib.pyplot as plt
import matplotlib as mpltlib
from matplotlib.colors import to_hex
from multiprocessing import Pool
from PAModelpy.PAModel import PAModel
import pandas as pd
import re
from scipy.optimize import minimize_scalar
from cobra import DictList
from collections import defaultdict
from typing import List, Dict, Literal, Optional
import warnings

from .PAM_data_classes import ValidationData, HyperParameters, ParametrizationResults, SectorConfig
from Modules.genetic_algorithm_parametrization import GAPOGaussian as GAPOGauss
from Modules.utils.error_calculation import calculate_r_squared_for_reaction, nanaverage
from Modules.utils.sampling_functions import adaptive_sampling
from Modules.utils.pam_generation import _extract_reaction_id_from_catalytic_reaction_id
from Modules.utils.sector_config_functions import (get_model_simulations_vs_sector,
                                                   perform_linear_regression,
                                                   change_sector_parameters_with_config_dict,
                                                   change_proteinsector_relation_from_growth_to_substrate_uptake,
                                                   SectorParameterDict)

#TYPEHINTS
Reaction2KcatDict = Dict[
    Literal['reaction', 'kcats', 'sensitivity'],
    Union[str, dict, float]
]
Enzyme2Reaction2KcatDict = Dict[str, List[Reaction2KcatDict]]

class PAMParametrizer():
    #from Schmidt et al(2016):
    TRANSLATIONAL_SECTOR_CONFIG = SectorConfig(
        sectorname='TranslationalProteinSector',
        slope= 0.04806975534478209, #g / g_cdw / h
        intercept=0.04738115630907698,#g / g_cdw
        substrate_range= list(np.arange(-4,0,1))
    )
    TRANSL_SECTOR_INTERCEPT_vs_GLC = 0.046136644909661115#g / g_cdw
    TRANSL_SECTOR_SLOPE_vs_GLC = -0.004340256958025938 #g / mmol_glc / h
    MEASURED_PROTEIN_FRACTION = 0.55*0.55
    MAXIMAL_INTERCEPT_UE = 0.67 #g_unusedprotein/g_p, as determined by O'Brien et al (2016)

    def __init__(self, pamodel:PAModel,
                 validation_data: Union[DictList[ValidationData], list, ValidationData],
                 hyperparameters: Optional[HyperParameters] = None,
                 sector_configs: Optional[Union[Dict[str, SectorConfig], bool]] = None,
                 substrate_uptake_id: str = "EX_glc__D_e",
                 max_substrate_uptake_rate: Union[float, int] = 0,
                 min_substrate_uptake_rate: Union[float, int] = -11,
                 sensitivity: bool = True,
                 minimal_unused_enzymes: float = 0.37,#g_ue/g_p
                 enzymes_to_evaluate:Optional[list] = None):

        self.core_genetic_algorithm = None

        #set up models with (find enzymes to change) and without (decrease simulation time) sensitivity
        self._pamodel = pamodel.copy(copy_with_pickle = True)
        self.pamodel_no_sensitivity = pamodel.copy(copy_with_pickle = True)
        self.pamodel_no_sensitivity.sensitivity = False

        #change total protein constraint to equality constraint for better fitting
        # self._set_total_protein_constraint_to_equality()
        if not hasattr(validation_data, "__iter__"): validation_data = [validation_data]
        self.validation_data = DictList(validation_data)
        self.hyperparameters = hyperparameters if hyperparameters is not None else HyperParameters()

        if sector_configs is None:
            sector_configs =  {
        'TranslationalProteinSector': self.TRANSLATIONAL_SECTOR_CONFIG}
        elif not sector_configs: #sector should not be changed
            sector_configs = {}

        self.sector_configs = sector_configs
        self.minimal_unused_enzymes = minimal_unused_enzymes

        self.substrate_uptake_ids = [csource_data.id for csource_data in validation_data]

        self.parametrization_results = ParametrizationResults(self.substrate_uptake_ids)
        self.parametrization_results.initiate_flux_results()
        self.enzyme_ids = [enzyme.id for enzyme in pamodel.enzymes]

        self.substrate_uptake_id = substrate_uptake_id
        self.min_substrate_uptake_rate = min_substrate_uptake_rate
        self.max_substrate_uptake_rate = max_substrate_uptake_rate
        self.fluxes_df = pd.DataFrame()
        self.final_error = 0
        self.previous_errors = 0
        self.error_fraction_to_deviate_between_runs = 0.1

        self.sensitivity = sensitivity
        self.enzymes_to_evaluate = enzymes_to_evaluate if enzymes_to_evaluate is not None else []

        # attributes for keeping track of the workflow
        self.iteration = 0
        self.bins = list(range(self.hyperparameters.number_of_bins))

        self._create_result_dirs()
        self.result_figure_file = os.path.join(
            "Results", "2_parametrization",
            "progress",
            f"pam_parametrizer_progress_{hyperparameters.filename_extension}.png")
        self.result_diagnostics_file = os.path.join(
            "Results",  "2_parametrization",
            "diagnostics",
            f"pam_parametrizer_diagnostics_{hyperparameters.filename_extension}.xlsx")
        self.calculate_sector_parameters_for_multiple_csources()

    @property
    def sector_parameters(self):
        return self.sector_configs

    @sector_parameters.setter
    def sector_parameters(self, sector_parameters: Dict[str, SectorConfig]):
        self.sector_configs = sector_parameters
        self.calculate_sector_parameters_for_multiple_csources(reset=True)

    def run(self, remove_subruns:bool = True, binned:str = "False") -> None:
        """ Run the parametrization framework.

        For each iteration of parametrization the following steps are taken:
        0. Initiate objects which are needed to store the results
        1. Optional: Split the range of substrate uptake rates to evaluate in a number of bins
        2. Optional: For each bin:
            i. Run simulations and calculate error to experimental measurements
            ii. Select the parameters to optimize
            iii. Run a genetic algorithm to optimize the selected parameters
        3. Select the enzymes to evaluate for the entire range of substrate uptake rates
        4. Initialize a genetic algorithm (Optional: use the results of the previous genetic algorithm runs)
        5. Run the genetic algorithm and use the best parameter set to reparametrize the model

        This is repeated until a max number of iterations is reached or until the simulation error is below a threshold

        Args:
            remove_subruns: If True, the results from the genetic algorithm runs within the bins are removed.
            The results of the final run won't be removed.
            binned: str: configuration of the binning process. There are 3 configurations possible:
                1. 'before': Only perform binning starting the first iteration of the workflow
                2. 'all': Perform binning for each iteration
                3. 'False'/other: No binning

        Returns:
            A png image of the parametrization progress and an excel file with the resulting parameters
        """
        self._set_pamodel_no_sensitivities()
        # self.calculate_sector_parameters_for_multiple_csources()

        #setup plot to visualize progress
        fig, axs = self.plot_valid_data()
        fig = self.plot_simulation(fig=fig, axs=axs, color="#010328", sensitivity = False)

        #keep track of time for computational performance
        start = time.perf_counter()
        if binned == "before":
            self._init_results_objects()
            files_to_remove = self.perform_iteration_in_bins(start)
            self.evaluate_and_save_results_of_iteration(start, files_to_remove, remove_subruns, fig, axs)

        #perform parametrization until a max number of iterations
        while (self.iteration+1 <= self.hyperparameters.threshold_iteration) & (self.final_error <= self.hyperparameters.threshold_error):
            # 0. initiate structures and objects to save results and keep track of process
            start_time_iteration = time.perf_counter()
            self.iteration += 1

            if binned == "all": #binning in all iterations
                # 1-4. Run simulations and genetic algorithm in bins and
                # restart a genetic algorithm using the resulting populations
                files_to_remove = self.perform_iteration_in_bins(start)
            else: # no binning
                # 1-4. Run genetic algorithm without bins
                files_to_remove = self.perform_iteration_without_bins()

            # 5. Reparametrize the model with best results and save/visualize results
            if files_to_remove is None:
                continue

            self.evaluate_and_save_results_of_iteration(start_time_iteration, files_to_remove, remove_subruns, fig, axs)

            if self._error_is_converging() and (self.iteration+1 <= self.hyperparameters.threshold_iteration):
                self.hyperparameters.number_of_kcats_to_mutate = len(self._pamodel.enzymes)
                print("Pick random enzymes to mutate because error is converging")
                self.iteration +=1
                files_to_remove = self.perform_iteration_without_bins()

                # files_to_remove = self.perform_iteration_in_bins(start)
                self.evaluate_and_save_results_of_iteration(start, files_to_remove, remove_subruns, fig, axs)


        self.optimize_sector_yintercept()
        self.save_final_diagnostics(figure = fig)
        plt.close(fig)

    def calculate_sector_parameters_for_multiple_csources(self, reset:bool = False):
        for vd in self.validation_data:
            for sector_id, sector_config in self.sector_configs.items():
                if sector_id in vd.sector_configs and not reset: continue
                sub_upt_id = vd.id

                #to prevent overflow metabolism, only select the lower growth rates to derive the equation
                substrate_range = self._get_substrate_range_lower_substrate_conc(vd.validation_range)
                vd.sector_configs[sector_id] = change_proteinsector_relation_from_growth_to_substrate_uptake(
                    pamodel = self.pamodel_no_sensitivity,
                    params = sector_config,
                    sector_id = sector_id,
                    substrate_uptake_id = sub_upt_id,
                    substrate_range = substrate_range
                )

    def perform_iteration_in_bins(self, start_time:float) -> list:
        self._init_results_objects()
        # 1. Get the binned substrate values and adjust binsize if required
        binned_substrates = self._bin_substrate_uptake_rates()


        # 2. Run model in bins, get sensitivities and calculate errors
        for subst_uptake_id, binned_substrate in binned_substrates.items():
            self._change_sector_parameters_for_substrate(substrate_uptake_id=subst_uptake_id,
                                                         pamodel=self._pamodel)
            for bin_id, bin_info in binned_substrate.items():
                self.process_bin(bin_id, bin_information=bin_info, substrate_uptake_id = subst_uptake_id)
                # print running time to check on progress
                print("time elapsed: ", time.perf_counter() - start_time, "sec, ", (time.perf_counter() - start_time) / 60, "min",
                      (time.perf_counter() - start_time) / 3600, "hour\n")

        #reset the validation dataframe for the entire range of substrate uptake rates
        self._init_validation_df()
        # 3. Restart genetic algorithm using populations from bins


        files_to_remove = self.restart_genetic_algorithm()

        return files_to_remove

    def perform_iteration_without_bins(self, random:bool=False) -> list:
        # 1. If there are no results yet, perform simulations for a range of glc_uptake_rates
        if self.iteration == 1:
            self._init_validation_df()
            self._init_results_objects()
            for substrate_uptake_rate in self.substrate_uptake_ids:
                self.run_simulations_to_plot(substrate_uptake_rate,
                                             save_fluxes_esc=True)

        # 2. Get the most sensitive enzymes from all runs
        if self.sensitivity and not random:
            enzymes_to_evaluate = self._determine_enzymes_to_evaluate_for_all_bins(
                nmbr_kcats_to_pick=self.hyperparameters.number_of_kcats_to_mutate)
        elif not self.sensitivity and not random:
            enzymes_to_evaluate = self.enzymes_to_evaluate
        else:
            enzymes_to_evaluate = self._get_random_enzymes_to_evaluate()


        # 3. run genetic algorithm for a range of substrate uptake rates
        self.run_genetic_algorithm(enzymes_to_evaluate=enzymes_to_evaluate,
                                   filename_extension=f"final_run_{self.iteration}")

        # 3. Restart genetic algorithm using populations from bins
        files_to_remove = self._get_genetic_algorithm_json_files()


        return files_to_remove

    def optimize_sector_yintercept(self,
                                  sector_id: str='UnusedEnzymeSector',
                                    throw_warning: bool = False
                                  ):
        """
            Optimize the y-intercept of a linear relation in a specified protein sector,
            keeping the x-intercept fixed. Updates sector parameters in validation_data.
            """
        with Pool() as pool:
            pool.starmap(self._optimize_sector_yintercept_for_validation_data,
                                       [(vd, sector_id, throw_warning) for vd in self.validation_data])


    def _optimize_sector_yintercept_for_validation_data(self,
                                                        vd: ValidationData,
                                                        sector_id: str,
                                                        throw_warning: bool):
        cache = {}
        def calculate_simulation_error(intercept):
            # avoid running multiple simulations with the same intercept
            rounded = round(intercept, 3)  # rounding to prevent float precision issues
            if rounded in cache:
                return cache[rounded]
            try:
                slope = -intercept / y0  # enforce x-intercept constraint
                self.pamodel_no_sensitivity.change_sector_parameters(sector,
                                                                     slope=slope,
                                                                     intercept=intercept,
                                                                     lin_rxn_id=vd.id
                                                                     )
                fluxes, substrate_rates = self.run_simulations_to_plot(vd.id,
                                                                       substrate_rates=vd.sampled_valid_data[
                                                                           f'{vd.id}_ub'],
                                                                       sensitivity=False)

                for substrate_rate, simulation_result in zip(substrate_rates, fluxes):
                    self.parametrization_results.add_fluxes_from_fluxdict(flux_dict=simulation_result,
                                                                          bin_id='intercept',
                                                                          substrate_reaction_id=vd.id,
                                                                          substrate_uptake_rate=substrate_rate,
                                                                          fluxes_abs=False)
                error = self._calculate_error_for_validation_data(valid_data=vd, bin_id='intercept')
                self.parametrization_results.remove_simulations_from_flux_df(substrate_uptake_id=vd.id,
                                                                             bin_id='intercept'
                                                                             )
                cache[rounded] = error
                return error
            except Exception as e:
                if throw_warning: warnings.warn(f"Simulation failed for intercept {intercept} on {vd.id}: {e}")
                return float('inf')  # penalize failure

        if not sector_id in vd.sector_configs: return

        sector_params = vd.sector_configs[sector_id].copy()
        y0 = sector_params['intercept'] / sector_params['slope']
        minimal_intercept = self.minimal_unused_enzymes * self._pamodel.total_protein_fraction
        res = minimize_scalar(
            calculate_simulation_error,
            bounds=(minimal_intercept, self.MAXIMAL_INTERCEPT_UE * self._pamodel.total_protein_fraction),
            # adjust upper bound as needed
            method='bounded',
            options={'xatol': 1e-3}
        )
        if not res.success:
            if throw_warning: warnings.warn(f"Optimization failed for {vd.id}: {res.message}")
            return
        sector_params['intercept'] = res.x
        sector_params['slope'] = res.x / y0

        print('new_sector parameters for', sector_id, vd.id)
        print(sector_params)
        vd.sector_configs[sector_id] = sector_params



    def evaluate_and_save_results_of_iteration(self, start_time_iteration:float, files_to_remove:list,
                                               remove_subruns:bool, fig: plt.Figure, axs: plt.Axes):

        self.reparametrize_pam()
        if self._pamodel_is_feasible:
            self._init_results_objects()
            # visualize results
            fig = self.plot_simulation(fig, axs, save_esc=True)

            for valid_data in self.validation_data:
                data = valid_data.valid_data
                substrate_uptake_id = valid_data.id
                substrate_rates_for_error = data[substrate_uptake_id+"_ub"].to_list()

                fluxes, substrate_rates = self.run_simulations_to_plot(substrate_uptake_id,
                                                                       substrate_rates=substrate_rates_for_error)

                for simulation_result, substrate_rate in zip(fluxes,substrate_rates):
                    self.parametrization_results.add_fluxes_from_fluxdict(flux_dict=simulation_result,
                                                                          bin_id="final",
                                                                          substrate_reaction_id= substrate_uptake_id,
                                                                          substrate_uptake_rate=substrate_rate,
                                                                          fluxes_abs=False)
            self.final_error = self.calculate_final_error()
            if pd.isna(self.final_error):
                print("Newly parametrized model was not feasible, reverting changes")
                self.revert_parametrization()
                return


            self.save_diagnostics(computational_time=time.perf_counter() - start_time_iteration)
            # 6. Display progress and repeat
            print("time elapsed: ", time.perf_counter() - start_time_iteration, "sec, ",
                  (time.perf_counter() - start_time_iteration) / 60, "min",
                  (time.perf_counter() - start_time_iteration) / 3600, "hour")
            print("Done with iteration number ", self.iteration)
            print("-------------------------------------------------------------------------------------------")

            if remove_subruns: self._remove_result_files(files_to_remove)
        else:
            print("Newly parametrized model is not feasible")
            return

    def process_bin(self,  bin_id: str, bin_information:list, substrate_uptake_id: str) -> None:
        """ Process a bin with specified substrate uptake rate information.

         Performs the following the following steps:
           2. Runs simulations in the range of substrate uptake rates using the protein allocation model.
           3. Calculates errors for different exchange rates.
           4. Checks for a solution; if none, the function returns early.
           5. Determines the most sensitive enzymes based on the calculated enzyme sensitivities.
           6. Determines if the bin should be split based on enzyme sensitivities.
           7. Runs the genetic algorithm on the bin to get mutated parameter value with a decreased simulation error

            Args:
                bin_information : A list containing substrate uptake rate information.
                                    Format: [start_rate, stop_rate, step_rate]
                bin_id : A unique identifier for the current bin.
           """

        print(
            f"The following range of substrate uptake rates will be analyzed: {bin_information[0]} - {bin_information[1]} for the {substrate_uptake_id} reaction")
        substrate_uptake_rates = self._init_validation_df(bin_information, [substrate_uptake_id])
        if substrate_uptake_rates is None: return
        for substrate_uptake_id, substrate_rates in substrate_uptake_rates.items():
            self.run_pamodel_simulations_in_bin(substrate_uptake_reaction=substrate_uptake_id,
                                                bin_id = bin_id,
                                                substrate_uptake_rates =substrate_rates
                                                )
            # calculate the error for the different exchange rates
            self.calculate_error(bin_id, substrate_uptake_id)
        # if there is no solution, continue to the next bin
        if len(self.parametrization_results.esc_df[self.parametrization_results.esc_df["bin"] == bin_id]) == 0:
            return
        # get the most sensitive enzymes
        top_enzyme_sensitivities = self.determine_most_sensitive_enzymes(bin_id,
                                                                         self.hyperparameters.number_of_kcats_to_mutate)
        self.determine_bin_to_split(top_enzyme_sensitivities, bin_id)
        enzymes_to_evaluate = self._parse_enzymes_to_evaluate(top_enzyme_sensitivities)
        self.run_genetic_algorithm(enzymes_to_evaluate =enzymes_to_evaluate,
                                   filename_extension=f"_iteration_{self.iteration}_bin_{bin_id}")

    def run_pamodel_simulations_in_bin(self,substrate_uptake_reaction:str,
                                       bin_id: Union[float, int], bin_information:list = None,
                                       substrate_uptake_rates: Union[list, np.array, pd.Series]=None) -> None:
        """ Run pamodel in range of substrate uptake rates defined in bin.

        Use the range of substrate uptake rate indicated in the bin_information to run simulations with the PAModel.
        Saves the simulation results in the parametrization_results data class.

        Args:
            substrate_uptake_reaction: id of the substrate uptake reactions for which to run the simulations
            bin_id: identifier of the bin to run
            bin_information: list with the substrate uptake range in the following order: [start, stop, step],
                    where start, stop and step relate to the start, end and stepsize of the substrate uptake rate range
            substrate_uptake_rates: iterable with the substrate uptake rates to test

        Notes:
            - Either bin_information or substrate_uptake_rates should be defined in order to define the substrate uptake
                    rates to constrain the model
        """
        self._change_sector_parameters_for_substrate(substrate_uptake_reaction, self._pamodel)
        if bin_information is not None:
            start, stop, step = bin_information[0], bin_information[1], bin_information[2]
            substrate_range = np.arange(start, stop, step)
        elif substrate_uptake_rates is not None:
            substrate_range = substrate_uptake_rates
            start = min(substrate_range)
            stop = max(substrate_range)
        print(f"The following range of substrate uptake rates will be analyzed: {start} - {stop}")
        for substrate_uptake_rate in substrate_range:
                self.run_pamodel_simulation(substrate_uptake_rate, bin_id, substrate_uptake_reaction)

    def run_pamodel_simulation(self, substrate_uptake_rate: Union[float, int],
                               bin_id: Union[str, float, int], substrate_uptake_reaction:str) -> None:
        """ Running pamodel simulation for specific substrate uptake rate.

        Running PAModel simulations and saving the resulting fluxes and enzymes sensitivities coefficient
        for later analysis

        Args:
            substrate_uptake_rate: used to constrain the PAModel
            bin_id: identifier of the bin in which this simulation is run (for saving purposes)
            substrate_uptake_reaction: id of the substrate uptake reaction
        """

        with self._pamodel:
            # change glucose uptake rate
            if substrate_uptake_rate < 0:
                self._pamodel.change_reaction_bounds(rxn_id=substrate_uptake_reaction,
                                                     lower_bound=substrate_uptake_rate, upper_bound=0)
            else:
                self._pamodel.change_reaction_bounds(rxn_id=substrate_uptake_reaction,
                                                     lower_bound=0, upper_bound=substrate_uptake_rate)
            # solve the model
            self._pamodel.optimize()
            if self._pamodel.solver.status == "optimal" and self._pamodel.objective.value != 0:
                self.save_pamodel_simulation_results(abs(substrate_uptake_rate), bin_id, substrate_uptake_reaction)

    def save_pamodel_simulation_results(self, substrate_uptake_rate: Union[float, int],
                                        bin_id: Union[str, float, int], substrate_uptake_reaction:str) -> None:
        """ Saving simulation results.

            Saving the resulting fluxes and enzymes sensitivities coefficient from a successful PAModel simulation
            for later analysis

            Args:
                substrate_uptake_rate: used to constrain the PAModel
                bin_id: identifier of the bin in which this simulation is run (for saving purposes)
                substrate_uptake_reaction: id of the substrate uptake reaction
        """

        flux_results = self.parametrization_results.flux_results.get_by_id(substrate_uptake_reaction)
        flux_results.substrate_range += [substrate_uptake_rate]

        self.parametrization_results.add_fluxes(self._pamodel, bin_id, substrate_uptake_reaction, substrate_uptake_rate,
                                                fluxes_abs=False)
        self.parametrization_results.add_enzyme_sensitivity_coefficients(self._pamodel.enzyme_sensitivity_coefficients,
                                                                         bin_id, substrate_uptake_rate)

    def calculate_error(self, bin_id: Union[float, int, str], substrate_uptake_reaction: str) -> None:
        """ Evaluate the model simulations compared to a reference dataset within a specified substrate uptake range.

        This method calculates the average difference between model simulations and reference data for the total substrate
        uptake range and available reactions within the specified bin.

        Notes:
        - The reference dataset is filtered to include only data points within the specified substrate uptake range.
        - Errors are calculated for different exchange rates and biomass reactions.
        - The 'calculate_r_squared_for_reaction' method is used to determine the R-squared value for each reaction.
        - The calculated errors are stored in the 'error_df' attribute for further analysis.

        Args:
            bin_id: Identifier for the bin under consideration.
            substrate_uptake_reaction: id of the substrate uptake reaction
        """
        valid_data = self.validation_data.get_by_id(substrate_uptake_reaction)
        reactions_to_validate = valid_data._reactions_to_validate
        validation_df = self._get_validation_data_to_validate(valid_data)
        # calculate error for different exchange rates
        error = self._calculate_error_for_reactions(substrate_uptake_id=substrate_uptake_reaction,
                                                    validation_df = validation_df,
                                                    reactions_to_validate = reactions_to_validate,
                                                    # substrate_range = self.validation_data.sampled_valid_data_df[self.substrate_uptake_id+'_ub'],
                                                    bin_id = bin_id)

        self.parametrization_results.add_error_to_error_df(substrate_uptake_reaction, bin_id, error)

    def calculate_final_error(self) -> float:
        """ Calculate the simulation error of the reparametrized model

        Makes use of the fluxes saved in self.parametrization_results to calculate the error

        Returns:
            final_error: mean of the r-squared value of the simulations for the different validation reactions
        """

        error = []
        for valid_data in self.validation_data:
            error += [self._calculate_error_for_validation_data(valid_data)]
            #remove the simulations from the result dataframe
            self.parametrization_results.remove_simulations_from_flux_df(valid_data.id,"final")
        final_error = np.nanmean(error)
        print("The final error is ", final_error)

        return final_error

    def _calculate_error_for_validation_data(self,
                                             valid_data:ValidationData,
                                             bin_id: str = 'final'
                                             ) -> float:
        substrate_uptake_id = valid_data.id
        reactions_to_validate = valid_data._reactions_to_validate
        validation_df = self._get_validation_data_to_validate(valid_data)

        return nanaverage(self._calculate_error_for_reactions(substrate_uptake_id=substrate_uptake_id,
                                                              validation_df=validation_df,
                                                              reactions_to_validate=reactions_to_validate,
                                                              bin_id=bin_id
                                                              ))

    def determine_most_sensitive_enzymes(self, bin_id: Union[str, float, int], nmbr_kcats_to_pick: int) -> pd.DataFrame:
        """
        Determine the top n (with n being the nmbr_of_kcats_to_pick) sensitive enzymes in a specific
        bin and save them to the dataclasses stored in parametrization_results. Selection is
        based on the enzyme sensitivity coefficients (refer to PAModelpy documentation for more
        details on calculations)

        Args:
            bin_id: bin identifier
            nmbr_kcats_to_pick: number of enzymes to select for optimization

        Returns:
            esc_topn_df: DataFrame containing the (mean) ESC values of the selected enzymes, with as columns:
                            [enzyme_id,  mean,  absolute_esc_mean,  bin,  substrate, rxn_id,  coefficient]
        """
        esc_results_df = self.parametrization_results.esc_df
        esc_in_bin = esc_results_df[esc_results_df["bin"] == bin_id]

        esc_topn_df = self._select_topn_enzymes(esc_results_df = esc_in_bin,
                                                nmbr_kcats_to_pick = nmbr_kcats_to_pick)
        return esc_topn_df

    def determine_bin_to_split(self, esc_topn_df: pd.DataFrame, bin_id: Union[str, float, int]) -> None:
        """ Determine whether a bin should be split based on the variability of ESC values.

        This function calculates the variability of ESC values within a bin and checks if it exceeds a certain threshold,
        indicating a significant change in the enzymatic substrate concentrations (ESCs) over the course of the bin.

        Args:
            esc_topn_df: DataFrame containing ESC values, where each column corresponds to an enzyme or substrate
            uptake rate, and rows represent different substrate uptake rates
            bin_id: Identifier for the bin under consideration.
        """
        enzyme_ids = esc_topn_df.enzyme_id.drop_duplicates()

        for enzyme_id in enzyme_ids:
            if enzyme_id == "substrate":
                continue

            esc_variability = self._calculate_esc_variability(esc_topn_df, enzyme_id)

            # Determine if ESC variability exceeds the threshold for splitting the bin
            if self._esc_variability_larger_than_threshold(esc_variability):
                self.parametrization_results.bins_to_change.loc[len(self.parametrization_results.bins_to_change)] = [bin_id, True, False]

    def run_genetic_algorithm(self,
                              enzymes_to_evaluate: Enzyme2Reaction2KcatDict,
                              filename_extension:str) -> None:
        """ Function to run the genetic algorithm.

        Results will be stored in a json, pickle and xlsx file with
        filename defined by: self.hyperparameters.genetic_algorithm_filename_base + filename_extension

        Args:
            enzymes_to_evaluate (dict): Dictionary with required information about which enzyme parameters should be optimized.
                Format:
                {
                    'enzyme_id': [{
                        'reaction': 'reaction_id',
                        'kcats': {'f':kcat_value_forward, 'b':kcat_value_reverse},
                        'sensitivity': esc_value
                    }]
                }
            filename_extension: extension of base filename to save the result of the genetic algorithm
        """
        substrate_rates = dict()
        sector_configs_per_substrate = dict()
        for valid_data in self.validation_data:
            #need to make sure that there is data sampled in the range to evaluate for this specific substrate
            if valid_data.sampled_valid_data is not None:
                substrate_rates[valid_data.id] = valid_data.sampled_valid_data[valid_data.id+"_ub"]
            sector_configs_per_substrate[valid_data.id] = valid_data.sector_configs

        ga = self._init_genetic_algorithm(substrate_uptake_rates=substrate_rates,
                                          enzymes_to_evaluate=enzymes_to_evaluate,
                                          sector_configs_per_substrate=sector_configs_per_substrate,
                                          filename_extension=filename_extension)
        ga.start()

    def restart_genetic_algorithm(self) -> None:
        json_files = self._get_genetic_algorithm_json_files()
        indiv_kcat_values_all_bins = pd.DataFrame()
        #Get the best individuals from the populations
        best_individuals = []
        for json in json_files:
            indiv_kcat_values_bin, _ = self._get_mutated_kcat_values_from_genetic_algorithm(json[:-4]+"xlsx")
            for i, row in indiv_kcat_values_bin.iterrows():
                best_individuals += {row.id:
                                         {"reaction": row.rxn_id, "direction": row.direction,
                                          "kcat": row.value, "sensitivity": 0.7}
                                     }
            if len(indiv_kcat_values_all_bins) == 0:
                indiv_kcat_values_all_bins = indiv_kcat_values_bin
            else:
                indiv_kcat_values_all_bins = pd.concat([indiv_kcat_values_all_bins, indiv_kcat_values_bin])

        enzymes_to_evaluate = self._determine_enzymes_to_evaluate_for_all_bins(
                                nmbr_kcats_to_pick= self.hyperparameters.number_of_kcats_to_mutate)
        if enzymes_to_evaluate is None: return None


        substrate_rates = {fr.id:fr.substrate_range for fr in self.parametrization_results.flux_results}
        sector_configs_per_substrate = {vd.id: vd.sector_configs for vd in self.validation_data}

        ga = self._init_genetic_algorithm(substrate_uptake_rates=substrate_rates,
                                          enzymes_to_evaluate=enzymes_to_evaluate,
                                          sector_configs_per_substrate=sector_configs_per_substrate,
                                          filename_extension=f"final_run_{self.iteration}")
        # TODO make a population with the current solution and the solution from the different ga's for the duplicate enzymes

        ga.restart(json_files)
        return json_files

    def reparametrize_pam(self, best_individual_kcat_df:pd.DataFrame=None)-> None:
        """ Change the selected kcat values to the optimized values

        Use the result of the genetic algorithm to change the parameters in the model if the error from the genetic
        algorithm is better than the error of the previous workflow iteration

        """
        if best_individual_kcat_df is None:
            best_individual_kcat_df, best_indiv_error = self._get_mutated_kcat_values_from_genetic_algorithm()

        # error_threshold = self.final_error*(1-self.error_fraction_to_deviate_between_runs)
        # if best_indiv_error < error_threshold and best_indiv_error>0:
        #     return
        for i, row in best_individual_kcat_df.iterrows():
            # current value is coefficient calculated as follows:
            # coeff = 1 / (kcat * 3600 * 1e-6) #3600 to convert to /h to /s *1e-6 to make calculations more accurate
            # need to convert the coefficient to the actual kcat
            new_kcat = 1 / (row["value"] * 3600 * 1e-6)
            direction = row["direction"]
            self._change_kcat_value_for_enzyme(enzyme_id= row["id"], kcat_dict = {row["rxn_id"]:
                                                                           {direction: new_kcat}
                                                                       })

    def revert_parametrization(self) -> None:
        for enzyme, rxn2kcat_info in self.enzymes_to_evaluate.items():
            for rxn2kcat in rxn2kcat_info:
                self._change_kcat_value_for_enzyme(enzyme_id=enzyme,
                                                   kcat_dict={
                                                       rxn2kcat['reaction']: {
                                                           d: 1/(kcat * 3600 * 1e-6) for d, kcat in rxn2kcat['kcats'].items()
                                                       }
                                                   })


    def save_diagnostics(self, computational_time: float,
                         results_filename:str = None) -> None:
        """
        Save the intermediate results to result objects

        Args:
            computational_time: time used for this iteration
            results_filename: filename where the results of the genetic algorithm are stored. If not giver, a default
                    filenmae is used
        """
        best_individual_kcat_df, best_indiv_error = self._get_mutated_kcat_values_from_genetic_algorithm(results_filename)

        for i, row in best_individual_kcat_df.iterrows():
            enz_rxn_kcat_info = self._parse_row_with_enz_rxn_kcat_for_saving(row)
            self.parametrization_results.add_best_individuals(run_id=self.iteration,
                                                              best_indiv_enz_rxn_kcat=enz_rxn_kcat_info,
                                                              ga_error= best_indiv_error)
        self.parametrization_results.add_computational_time(run_id=self.iteration,
                                                            time_in_sec= computational_time)

        self.parametrization_results.add_final_error(run_id= self.iteration,
                                                     final_error= self.final_error)

    def save_final_diagnostics(self, figure:plt.Figure = None) -> None:
        """Save the results of all iterations to excel and png

        Args:
            figure: figure to save (progress figure)
        """
        if figure is not None:
            figure.savefig(self.result_figure_file, dpi=100, bbox_inches="tight")
        sector_configs = pd.DataFrame(columns = ["substrate_uptake_id", "slope", "intercept", 'sector_id'])
        ga_weights = pd.DataFrame({'reaction':self.hyperparameters.genetic_algorithm_hyperparams['error_weights'].keys(),
                                   'weight': self.hyperparameters.genetic_algorithm_hyperparams['error_weights'].values()} )

        for vd in self.validation_data:
            for sector_id, sector_config in vd.sector_configs.items():
                sector_configs = pd.concat(
                    [sector_configs, pd.Series(
                        {
                            "substrate_uptake_id":vd.id,
                            "slope": sector_config['slope'],
                            "intercept": sector_config['intercept'],
                            "sector_id": sector_id
                        }).to_frame().T],
                    ignore_index=True
                )


        with pd.ExcelWriter(self.result_diagnostics_file) as writer:
            # Write each DataFrame to a specific sheet
            self.parametrization_results.best_individuals.to_excel(writer, sheet_name="Best_Individuals", index=False)
            self.parametrization_results.computational_time.to_excel(writer, sheet_name="Computational_Time",
                                                                     index=False)
            self.parametrization_results.final_errors.to_excel(writer, sheet_name="Final_Errors", index=False)
            sector_configs.to_excel(writer, sheet_name="sector_parameters", index = False)
            if len(ga_weights)>0:ga_weights.to_excel(writer, sheet_name="reaction_weights", index = False)

    def add_new_substrate_source(self, new_substrate_uptake_id:str,
                                 validation_data: pd.DataFrame,
                                 substrate_range: list[Union[int, float]],
                                 reactions_to_validate: list):
        valid_data = ValidationData(validation_data,
                                    new_substrate_uptake_id,
                                    substrate_range)
        valid_data._reactions_to_validate = reactions_to_validate
        # if new_substrate_uptake_id not in self.substrate_uptake_ids:
        #     self.substrate_uptake_ids.append(new_substrate_uptake_id)
        self.validation_data.append(valid_data)
        self.parametrization_results.add_new_substrate_source(new_substrate_uptake_id,
                                                              reactions_to_validate)
    @property
    def pamodel(self) -> PAModel:
        return self._pamodel

    @pamodel.setter
    def pamodel(self, pamodel:PAModel):
        self._pamodel = pamodel
        self._set_pamodel_no_sensitivities()

    ###########################################################################################################
    #WORKER FUNCTIONS
    ###########################################################################################################
    def _set_pamodel_no_sensitivities(self) -> None:
        self.pamodel_no_sensitivity = self._pamodel.copy(copy_with_pickle=True)
        self.pamodel_no_sensitivity.sensitivity = False

    def _create_result_dirs(self, result_file_path: str = os.path.join(
            "Results", "2_parametrization")) -> None:
        """ Check if directories for storing the results exist.

        If the directories do not exist, they are created

        Args:
            result_file_path: the path to the directory where the results should be stored
        """
        if not os.path.isdir(result_file_path):
            os.mkdir(result_file_path)
        for subdir in ['progress', 'diagnostics']:
            if not subdir in os.listdir(result_file_path):
                os.mkdir(os.path.join(result_file_path, subdir))

    def _change_sector_parameters_for_substrate(self,
                                                substrate_uptake_id: str,
                                                pamodel:PAModel
                                                ) -> None:
        sector_configs = self.validation_data.get_by_id(substrate_uptake_id).sector_configs
        for sector_id, sector_config in sector_configs.items():
            change_sector_parameters_with_config_dict(pamodel,
                                                      sector_config=sector_config,
                                                      substrate_uptake_id = substrate_uptake_id,
                                                      sector_id = sector_id
                                                      )

    def _get_substrate_range_lower_substrate_conc(self,
                                                  validation_range:list[Union[int, float]],
                                                  number_of_steps: int = 5
                                                  ) -> Iterable:
        #only loop over the low growth rates, to prevent overflow like metabolism to interfere with the derivation of a linear equation
        # Determine the range to loop over based on absolute values
        start = min(abs(validation_range[0]), abs(validation_range[1])) / 2
        end = max(abs(validation_range[0]), abs(validation_range[1])) / 2
        step = abs(start-end)/number_of_steps

        # Create the substrate range
        substrate_range = np.arange(start, end, step)

        # Adjust the sign if the original range was negative
        if validation_range[0] <= 0 and validation_range[1] <= 0:
            substrate_range = -substrate_range[::-1]
        return substrate_range

    def _pamodel_is_feasible(self):
        #check for feasibility at a mid substrate uptake rate
        model = self.pamodel_no_sensitivity
        sub_upt = self.max_substrate_uptake_rate -self.min_substrate_uptake_rate/2
        if sub_upt<0:
            model.change_reaction_bounds(self.substrate_uptake_id, lower_bound=sub_upt, upper_bound=0)
        else:
            model.change_reaction_bounds(self.substrate_uptake_id, lower_bound=0, upper_bound=sub_upt)
        model.optimize()
        return model.solver.status == "optimal"

    def _bin_substrate_uptake_rates(self):
        """ Bin the substrate uptake rate in intervals.

         If required, the binsize is adjusted

        Returns:
            dict: dictionary with bin_id: [start, end, stepsize] key value pairs
        """
        binned_substrates = {}
        for valid_data in self.validation_data:
            substrate_start = valid_data.validation_range[0]
            binned_substrate = {}
            bin_range = abs((abs(valid_data.validation_range[0]) - abs(valid_data.validation_range[1])))/self.hyperparameters.number_of_bins
            for i in self.bins:
                new_bin = self._make_new_bin(bin_id = i, bin_range = bin_range, substrate_start = substrate_start)
                binned_substrate = {**binned_substrate, **new_bin}
                #update starting concentration for new bin
                substrate_start += bin_range
            binned_substrates[valid_data.id]= binned_substrate
        return binned_substrates

    def _make_new_bin(self, bin_id:Union[float, int], bin_range: Union[float, int],
                      substrate_start:Union[float, int]):
        """ Creating a dictionary with all information about a single bin (range of substrate uptake rates).

        Saves the start, end and stepsize per bin.
        The function also checks if the bin should be adapted based on the previous iteration of the workflow

        Args:
            bin_id: identifier of the bin to make
            bin_range: length of the bin in the units of the substrate uptake rate
            substrate_start: the first substrate uptake rate in the bin

        Returns:
            dict: dictionary with bin_id: [start, end, stepsize] key value pairs
        """
        # if bins should be smaller, adjust binsize
        if ((len(self.parametrization_results.bins_to_change) > 0) and
                (bin_id in self.parametrization_results.bins_to_change["bin"].values)):
            new_bin = self._adjust_binsize(bin_id=bin_id, bin_range=bin_range, substrate_start=substrate_start)
        else:
            stepsize = bin_range / self.hyperparameters.bin_resolution
            new_bin = {bin_id: [substrate_start, substrate_start + bin_range, stepsize]}
        return new_bin

    def _adjust_binsize(self, bin_id, bin_range, substrate_start):
        bins_to_change = self.parametrization_results.bins_to_change
        to_split = [(bins_to_change[bins_to_change["bin"] == bin_id]) & (bins_to_change["split"])]
        if len(to_split) > 0:
            stepsize = bin_range * 0.5 / self.hyperparameters.bin_resolution
            return {**{bin_id: [substrate_start, substrate_start + bin_range * 0.5, stepsize]},
                                **{bin_id + 0.1: [substrate_start + bin_range * 0.5, substrate_start + bin_range, stepsize]}}


    def _init_genetic_algorithm(self,
                                substrate_uptake_rates: dict,
                                enzymes_to_evaluate: Enzyme2Reaction2KcatDict,
                                sector_configs_per_substrate: Dict[str, SectorParameterDict],
                                filename_extension: str
                                ) -> GAPOGauss:
        """
        Initializes the core genetic algorithm object.

        Upon initiation, the core genetic algorithm initializes the genetic algorithm runner (ga_param) and the fitness
        evaluation class. Hyperparameters related to the genetic algorithm can be adjusted in
        `self.hyperparameters.genetic_algorithm_hyperparameters` (dict). The sampling method (uniform or from a normal
        distribution) can be chosen by changing the core genetic algorithm object in `self.hyperparameters.genetic_algorithm`.

        Files will be saved with a filename defined by:
        `self.hyperparameters.genetic_algorithm_filename_base + filename_extension`

        Args:
            substrate_uptake_rates (list): Information about the substrate uptake rates to evaluate in the form of [start,..,stop].
            enzymes_to_evaluate (dict): Dictionary with required information about which enzyme parameters should be optimized.
                Format:
                {
                    'enzyme_id': [{
                        'reaction': 'reaction_id',
                        'kcats': {'f':kcat_value_forward, 'b':kcat_value_reverse},
                        'sensitivity': esc_value
                    }]
                }
            sector_configs_per_substrate (dict of dict): Dictionary with the slope and intercept of the translational and unused
                sector configuration for each substrate.
                Format:
                {'substrate_uptake_id':{'ProteinSector'{
                    'slope':float, #slope in g/mmol/h
                    'intercept':float #intercept in g/mmol
                    }}
                }
            filename_extension (str): Extension of the base filename to save the result of the genetic algorithm.

        Returns:
            ga: Core genetic algorithm object ready to be started.
        """

        results_filename = self.hyperparameters.genetic_algorithm_filename_base + filename_extension

        ga = self.hyperparameters.genetic_algorithm(
            model=self.pamodel_no_sensitivity.copy(copy_with_pickle=True),
            enzymes_to_eval= enzymes_to_evaluate,
            sector_configs_per_substrate = sector_configs_per_substrate,
            valid_data = self._create_validation_data_dict_for_genetic_algorithm(),#{valid_data.id: valid_data.sampled_valid_data[valid_data._reactions_to_validate + [valid_data.id+"_ub"]] for valid_data in self.validation_data if valid_data.sampled_valid_data is not None},
            filename_save = results_filename,
            substrate_uptake_id = self.substrate_uptake_id,
            substrate_uptake_rates = substrate_uptake_rates, # bin_info: [start, stop, step]
            objective_id = self.validation_data[0]._get_biomass_reactions(),
            **self.hyperparameters.genetic_algorithm_hyperparams
        )
        return ga

    def _create_validation_data_dict_for_genetic_algorithm(self) -> Union[dict, pd.DataFrame]:
        valid_data_dict = {}
        for valid_data in self.validation_data:
            df_for_validation = self._get_validation_data_to_validate(valid_data)
            if df_for_validation is not None:
                valid_data_dict[valid_data.id] = df_for_validation
        return valid_data_dict

    def _get_validation_data_to_validate(self, valid_data:ValidationData) -> pd.DataFrame:
        if valid_data.sampled_valid_data is not None:
            columns_to_sample = ([rxn for rxn in valid_data._reactions_to_validate
                                 if rxn in valid_data.sampled_valid_data]  #need to make sure there is data available
                                 + [valid_data.id + "_ub"]) #also check the predicted substrate uptake rate
            return valid_data.sampled_valid_data[columns_to_sample]

    def _init_validation_df(self, bin_information:list = None, substrate_uptake_ids: list = None) -> dict:
        # Only sample those points which are lower than upper bound and higher than lower bound
        # Also correct for reaction direction
        substrate_rates = {}
        if bin_information is not None:
            # if the range is provided, make sure to get the right directionality based on whether the sign is positive or negative
            lower, upper = self._correct_upper_and_lower_ranges_validation_data_for_sign(bin_information)

        validation_data_objects = self.validation_data
        if substrate_uptake_ids is not None:
            validation_data_objects = [validation_data_objects.get_by_id(sub_id) for sub_id in substrate_uptake_ids]

        for validation_data in validation_data_objects:
            # get the ranges to validate for the individual carbon sources if the range is not provided
            if bin_information is None:
                lower, upper = self._correct_upper_and_lower_ranges_validation_data_for_sign(
                    validation_data.validation_range)
            substrate_uptake_id = validation_data.id

            validation_df = validation_data.valid_data[
                 (abs(validation_data.valid_data[substrate_uptake_id + "_ub"]) >= lower) &
                 (abs(validation_data.valid_data[substrate_uptake_id + "_ub"]) <= upper)
             ]
            #no data to validate within this range available

            if len(validation_df) == 0: continue
            validation_data.sampled_valid_data = adaptive_sampling(exp_data = validation_df,
                                                            num_samples=self.hyperparameters.bin_resolution)
            validation_data.sampled_valid_data.reset_index()
            substrate_rates[substrate_uptake_id] = validation_df[substrate_uptake_id + "_ub"]

        #no substrate uptake rates in the requested range
        if len([rate for rate in substrate_rates.values()]) == 0:return None
        return substrate_rates
        # get reference datapoints

    def _correct_upper_and_lower_ranges_validation_data_for_sign(self, validation_range:list):
        if validation_range[0] < 0:
            upper = abs(validation_range[0])
            lower = abs(validation_range[1])
        else:
            upper = validation_range[1]
            lower = validation_range[0]
        return lower, upper


    def _init_results_objects(self):
        reactions_to_validate_dict = {}
        for rxn in self.parametrization_results.substrate_uptake_reactions:
            reactions_to_validate_dict[rxn] = self.validation_data.get_by_id(rxn)._reactions_to_validate

        self.parametrization_results.initiate_result_dfs(reactions_to_validate=reactions_to_validate_dict)


    def _select_topn_enzymes(self, esc_results_df: pd.DataFrame, nmbr_kcats_to_pick:int) -> pd.DataFrame:
        """
                Determines the top n (with n being `nmbr_of_kcats_to_pick`) sensitive enzymes.

                Selection is based on the enzyme sensitivity coefficients (refer to PAModelpy documentation for more details on calculations).

                Args:
                    esc_results_df (DataFrame): Dataframe with enzyme sensitivity coefficients obtained from the PAModel simulations.
                    nmbr_kcats_to_pick (int): Number of enzymes to select for optimization.

                Returns:
                    esc_topn_df (DataFrame): DataFrame containing the (mean) ESC values of the selected enzymes, with columns:
                        [enzyme_id, mean, absolute_esc_mean, bin, substrate, rxn_id, coefficient]
                """
        # 0. get data in shape: make sure each reaction has their own row and substrate value is absolute
        # Splitting rxn_id column and creating new rows
        esc_results_df["rxn_id"] = esc_results_df["rxn_id"].str.split(",")
        esc_results_df = esc_results_df.explode("rxn_id", ignore_index=True)
        esc_results_df["substrate"] = esc_results_df["substrate"].abs()

        #make sure each gpr is only present once
        esc_results_df['reaction_parsed'] = esc_results_df.rxn_id.apply(
            _extract_reaction_id_from_catalytic_reaction_id
        )
        esc_results_df = esc_results_df.drop_duplicates(
            ['enzyme_id', 'reaction_parsed']
        ).drop('reaction_parsed', axis = 1)

        # 1. determine top n values based on average ESC
        # Group by enzyme and calculate the average coefficient for each enzyme
        esc_results_grouped = esc_results_df.groupby("enzyme_id")["coefficient"].agg(["mean"]).reset_index()
        esc_results_grouped["absolute_esc_mean"] = esc_results_grouped["mean"].abs()
        # Sort the data based on the average coefficient in descending order
        esc_results_sorted = esc_results_grouped.sort_values(by="absolute_esc_mean", ascending=False)
        # Select the top n enzymes
        top_n_enzymes = esc_results_sorted.head(nmbr_kcats_to_pick)

        # 2. merge with the original DataFrame to connect to the information
        esc_topn_df = pd.merge(top_n_enzymes, esc_results_df.drop_duplicates(), on="enzyme_id")
        # get ids and select those from the dataframe
        esc_topn_df = esc_topn_df.sort_values(by="substrate")
        return esc_topn_df

    def _calculate_esc_variability(self, esc_df:pd.DataFrame ,enzyme_id: str) -> float:
        """
        Calculates ESC variability as a fractional change.

        Args:
            esc_df (DataFrame): DataFrame containing ESC values.
            enzyme_id (str): ID of the enzyme.

        Returns:
            float: ESC variability as a percentage change.
        """

        start_esc = esc_df[esc_df["enzyme_id"] == enzyme_id].coefficient.iloc[0]
        end_esc = esc_df[esc_df["enzyme_id"] == enzyme_id].coefficient.iloc[-1]

        # check if the ESCs are nonzero
        if start_esc != 0:
            esc_variability = abs((end_esc - start_esc) / start_esc)
        else:
            # Handle the case when start_esc is zero to avoid division by zero
            esc_variability = 0 if end_esc == 0 else 1  # Assuming that if start_esc is zero, any change is considered 100% change.

        return esc_variability

    def _esc_variability_larger_than_threshold(self, esc_variability: float)-> bool:
        return (esc_variability >= self.hyperparameters.bin_split_deviation_threshold)

    def _get_random_enzymes_to_evaluate(self):
        if self._pamodel_is_feasible():
            esc_df = self._pamodel.enzyme_sensitivity_coefficients
            esc_df["rxn_id"] = esc_df["rxn_id"].str.split(",")
            esc_df = esc_df.explode("rxn_id", ignore_index=True)
            enzymes_to_evaluate = self._parse_enzymes_to_evaluate(
                esc_df.sample(
                    self.hyperparameters.number_of_kcats_to_mutate
                ).rename({"coefficient": "mean"}, axis=1)
            )
        else:
            fake_esc_df = pd.DataFrame(columns=['enzyme_id', 'mean', 'absolute_esc_mean', 'bin', 'substrate', 'rxn_id', 'coefficient'])
            all_proteins = self._pamodel.enzyme_variables.copy()
            sampled_proteins = [self._pamodel.enzymes.get_by_id(
                ev.id) for ev in random.sample(all_proteins, k=self.hyperparameters.number_of_kcats_to_mutate)]
            for protein in sampled_proteins:
                rxn = random.sample(list(protein.rxn2kcat.keys()), k=1)[0]
                fake_esc_df.loc[len(fake_esc_df)] = [protein.id, 0,0,'',0,rxn,0]
            enzymes_to_evaluate = self._parse_enzymes_to_evaluate(
                fake_esc_df.sample(
                    self.hyperparameters.number_of_kcats_to_mutate
                ).rename({"coefficient": "mean"}, axis=1)
            )
        return enzymes_to_evaluate

    def _calculate_error_for_reactions(self, substrate_uptake_id:str,
                                       validation_df: pd.DataFrame,
                                       reactions_to_validate: list,
                                        bin_id: Union[float, int] = None) -> float:
        # calculate error for different exchange rates
        error = []
        flux_df = self.parametrization_results.flux_results.get_by_id(substrate_uptake_id).fluxes_df

        if len(flux_df) == 0:  # means model is infeasible
            return np.nan
        # check if we want to calculate the error for a single bin
        if bin_id is not None: flux_df = flux_df[flux_df["bin"] == bin_id]

        for rxn in reactions_to_validate:
            # only select the rows which are filled with data
            if rxn not in validation_df.columns: continue
            validation_data = validation_df.dropna(axis=0, subset=rxn)
            # if there are no reference data points, continue to the next reaction
            if len(validation_data) == 0:
                error += [np.nan]
                continue

            r_squared = calculate_r_squared_for_reaction(rxn, validation_df, substrate_uptake_id,
                                                               flux_df)

            error += [r_squared]
        return error


    def _parse_enzymes_to_evaluate(self, esc_topn_df: pd.DataFrame) -> Reaction2KcatDict:
        """
        Parses the enzyme-reaction pairs for which the kinetic parameters should be adjusted from the enzyme sensitivity coefficients.

        Args:
            esc_topn_df (DataFrame): DataFrame containing the (mean) ESC values of the selected enzymes, with columns:
                [enzyme_id, mean, absolute_esc_mean, bin, substrate, rxn_id, coefficient]

        Returns:
            dict: Dictionary with required information about which enzyme parameters should be optimized. Format:
                {
                    'enzyme_id': {
                        'reaction': 'reaction_id',
                        'kcat': kcat_value,
                        'sensitivity': esc_value
                    }
                }
        """
        enzymes_to_evaluate = defaultdict(list)
        for index, row in esc_topn_df.iterrows():
            enzyme_id = row["enzyme_id"]
            rxn_id = row["rxn_id"]
            enzyme_dict = {}
            enzyme_dict["sensitivity"] = row["mean"]
            #create the connection to the catalytic reaction related to the enzyme
            rxn_id = _extract_reaction_id_from_catalytic_reaction_id(rxn_id) #make sure the actual reaction id is used
            rxn_id = f"CE_{rxn_id}_{enzyme_id}"

            enzyme_dict["reaction"] = rxn_id
            kcat_dict = self._pamodel.enzymes.get_by_id(enzyme_id).rxn2kcat[rxn_id]
            # kcats are in 1/s needs to be converted to 1/h.
            # In the genetic algorithm, the kcats are directly set as coefficients
            # (1/(kcat *1e-6), where 1e-6 is used to scale enzyme concentrations)
            enzyme_dict["kcats"] = {dir: 1/(kcat* 3600 * 1e-6) for dir, kcat in kcat_dict.items() if kcat>0}
            enzymes_to_evaluate[enzyme_id].append(enzyme_dict)
        return dict(enzymes_to_evaluate)

    def _get_genetic_algorithm_json_files(self, subset:str = ""):
        filename_base = self.hyperparameters.genetic_algorithm_filename_base
        folder_path = self.hyperparameters.genetic_algorithm_hyperparams["folderpath_save"]

        # Define the regular expression pattern for matching file names
        file_pattern = re.compile(fr"{filename_base}{subset}.*\.json$")
        # Get a list of all files in the folder and filter using the regular expression
        json_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if file_pattern.match(f)]
        return json_files

    def _determine_enzymes_to_evaluate_for_all_bins(self, nmbr_kcats_to_pick) -> dict:
        """
        Parses enzyme-reaction pairs for adjusting kinetic parameters from enzyme sensitivity coefficients.

        Args:
            nmbr_kcats_to_pick (int): Number of enzymes to select based on the enzyme sensitivity coefficients.

        Returns:
            dict: Dictionary with required information about which enzyme parameters should be optimized. Format:
                {
                    'enzyme_id': {
                        'reaction': 'reaction_id',
                        'kcat': kcat_value,
                        'sensitivity': esc_value
                    }
                }
        """
        esc_results_df = self.parametrization_results.esc_df
        if len(esc_results_df) == 0: return self._get_random_enzymes_to_evaluate() #if infeasible get default parameter set to adjust
        esc_topn_df = self._select_topn_enzymes(esc_results_df,
                                                nmbr_kcats_to_pick)
        enzymes_to_evaluate = self._parse_enzymes_to_evaluate(esc_topn_df)
        self.enzymes_to_evaluate = enzymes_to_evaluate
        return enzymes_to_evaluate

    def _get_mutated_kcat_values_from_genetic_algorithm(self, results_filename:str = None):
        if results_filename is None:
            filename_base = self.hyperparameters.genetic_algorithm_filename_base
            folder_path = self.hyperparameters.genetic_algorithm_hyperparams["folderpath_save"]

            results_filename = os.path.join(folder_path,
                                                 filename_base + f"final_run_{self.iteration}.xlsx")
        final_run_results_best_indiv = pd.read_excel(results_filename, sheet_name="best_individual").drop(["type"],
                                                                                                          axis = 1)
        error_df = pd.read_excel(results_filename, sheet_name= "final_population")
        error = float(error_df.iloc[0]["fitness_weighted_sum"])
        return final_run_results_best_indiv ,error

    def _change_kcat_value_for_enzyme(self, enzyme_id:str, kcat_dict:dict) -> None:
        for pamodel in [self._pamodel, self.pamodel_no_sensitivity]:
            pamodel.change_kcat_value(enzyme_id=enzyme_id, kcats=kcat_dict)

    def _remove_result_files(self, file_base: Union[list, str]) -> None:
        """ Removes files resulting from genetic algorithm runs.

        Args:
            file_base (list, string): list of stings or pathlib.Path objects or a string with the files
            to remove. If there are file extensions, they will be removed.
        """
        if not hasattr(file_base, "__iter__"): file_base = [file_base]
        for file in file_base:
            if not isinstance(file, str): file = str(file)
            file_path_base = file.split(".")[0]
            [os.remove(file_path_base + file_type) for file_type in [".json", ".xlsx", ".pickle"]]

    def _set_total_protein_constraint_to_equality(self):
        for pamodel in [self._pamodel, self.pamodel_no_sensitivity]:
            pamodel.constraints[pamodel.TOTAL_PROTEIN_CONSTRAINT_ID].lb = pamodel.constraints[pamodel.TOTAL_PROTEIN_CONSTRAINT_ID].ub

    def _parse_row_with_enz_rxn_kcat_for_saving(self, enz_rxn_kcat_row:pd.Series) -> list:
        """ Parsing a row from a result dataframe for saving

        The row should contain the reaction id, enzyme id and kcat value.

        Args:
            enz_rxn_kcat_row: pandas.Series with 'id' (enzyme id), 'rxn_id' (reaction id) and
            'value' (kcat value in h-1) as keys.

        Returns:
            result: list with the following structure: [reaction_id, enzyme_id, kcat_in_s-1]

        """
        #need to convert the units of the kcat from 1/h to 1/s for easy interpretation
        enz_rxn_kcat_row["value"] = 1/(enz_rxn_kcat_row["value"]*3600*1e-6)
        return enz_rxn_kcat_row.to_list()

    def _error_is_converging(self):
        errors = self.parametrization_results.final_errors.copy()
        errors["diff"] = errors["r_squared"].diff().abs()
        return errors["diff"][-self.hyperparameters.threshold_nmbr_convergence:].max() < self.hyperparameters.threshold_convergence


    #################################################################################################################
    # VISUALIZATION
    #################################################################################################################

    def plot_valid_data(self):
        # plot flux changes with glucose uptake
        reactions_to_plot = self.validation_data.get_by_id(self.substrate_uptake_id)._reactions_to_plot

        num_reactions = len(reactions_to_plot) + 1
        num_cols = int(np.ceil(np.sqrt(num_reactions)))
        num_rows = int(np.ceil(num_reactions / num_cols))

        fig, axs = plt.subplots(num_rows,num_cols, dpi=100)

        for r, ax in zip(reactions_to_plot, axs.flatten()[:-1]):
            valid_data = self.validation_data.get_by_id(self.substrate_uptake_id)
            # plot data
            x = [abs(sub) for sub in valid_data.valid_data[self.substrate_uptake_id +"_ub"]]
            y = [abs(data) for data in valid_data.valid_data[r]]
            ax.set_ylabel(r)
            ax.scatter(x, y,
                           color="black", marker="o", s=30, linewidths=1.3,
                           facecolors=None, zorder=0,
                           label="Data")
            ax.set_xlabel(self.substrate_uptake_id + " $mmol/g_{CDW}/h$")

        # Remove any unused subplots
        for j in range(len(reactions_to_plot), len(axs)):
            fig.delaxes(axs[j])

        #make the plots square to easily see the relation
        axs.flatten()[-1].set_aspect("equal")
        axs.flatten()[-1].set_ylabel("Experimental measured fluxes")
        axs.flatten()[-1].set_xlabel("simulated fluxes")

        plt.ion()   # set interactive mode
        fig.tight_layout()
        fig.show()

        return fig, axs

    def plot_simulation(self, fig, axs,
                        return_fluxes:bool = False,
                        save_esc = False,
                        color:int= None,
                        cbar_label:str = "Iteration",
                        sensitivity = True) -> plt.Figure:

        if color is None:
            #adjust color to visualize progress
            #get viridis color palette
            cmap = plt.get_cmap("viridis")
            color = to_hex(cmap(self.iteration / self.hyperparameters.threshold_iteration))
            norm = mpltlib.colors.BoundaryNorm(list(range(self.hyperparameters.threshold_iteration+1)), cmap.N)

        fluxes_dict = {}
        substrate_range_dict = {}

        for substrate_id in self.substrate_uptake_ids:
            if substrate_id == self.substrate_uptake_id: substrate_range = None
            else: substrate_range = self.validation_data.get_by_id(substrate_id).valid_data[f'{substrate_id}_ub']


            fluxes, substrate_range = self.run_simulations_to_plot(substrate_uptake_id=substrate_id,
                                                                   substrate_rates=substrate_range,
                                                                    save_fluxes_esc=save_esc,
                                                                   sensitivity = sensitivity)
            if len(fluxes) > 0: # only plot feasible model results
                fluxes_dict[substrate_id] = fluxes
                substrate_range_dict[substrate_id] = substrate_range


        alpha = 1
        #only plot detailed info for the 'main' substrate
        if len(substrate_range_dict)!=0: #only plot if the model is feasible
            for r, ax in zip(self.validation_data.get_by_id(self.substrate_uptake_id)._reactions_to_plot, axs.flatten()[:-1]):
                # plot data
                line = ax.plot(substrate_range_dict[self.substrate_uptake_id].values(),
                               [abs(f[r]) for f in fluxes_dict[self.substrate_uptake_id]], linewidth=2.5,
                               zorder=5, color=color, alpha=alpha)
        # Plot a flux comparison to experimental data for the other fluxes
        for substr_id, fluxes in fluxes_dict.items():
            valid_data = self.validation_data.get_by_id(substr_id)
            if substr_id != self.substrate_uptake_id and valid_data.sampled_valid_data is not None:
                feas_substrate_rates = substrate_range_dict[substr_id].keys()
                feas_sampled_data = valid_data.valid_data#[
                    #valid_data.valid_data[substr_id + "_ub"].isin(feas_substrate_rates)
                #]
                for reaction in valid_data._reactions_to_plot:
                    exp_measurements = feas_sampled_data[reaction]
                    simulations = [f[reaction] for f in fluxes]
                    axs.flatten()[-1].scatter(exp_measurements, simulations, color = color, alpha = alpha)

        # Add colorbar
        if self.iteration == 1:
            cbar = fig.colorbar(mpltlib.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axs)
            cbar.set_label(cbar_label)

        fig.canvas.draw()
        fig.canvas.flush_events()
        if return_fluxes: return fig, fluxes_dict, substrate_range_dict
        return fig

    def run_simulations_to_plot(self, substrate_uptake_id:str,
                                substrate_rates: Union[np.array, list, pd.Series] = None,
                                save_fluxes_esc:bool = False,
                                sensitivity: bool = True) -> Tuple[list, dict[float,float]]:
        fluxes = list()
        substrate_range = dict()

        if sensitivity:
            pamodel = self._pamodel
        else:
            pamodel = self.pamodel_no_sensitivity

        previous_bounds = pamodel.get_reaction_bounds(substrate_uptake_id)

        #change the sector parameters for correct relation between substrate uptake rate, ribosome content and unused proteins
        self._change_sector_parameters_for_substrate(substrate_uptake_id,
                                                     pamodel)

        if substrate_rates is None:
            step = (self.max_substrate_uptake_rate-self.min_substrate_uptake_rate)/10
            substrate_rates = np.arange(self.min_substrate_uptake_rate, self.max_substrate_uptake_rate, step)
        for substrate in substrate_rates:
            if substrate>=0:
                pamodel.change_reaction_bounds(rxn_id=substrate_uptake_id,
                                                lower_bound=0, upper_bound=substrate)
            else:
                pamodel.change_reaction_bounds(rxn_id=substrate_uptake_id,
                                                    lower_bound=substrate, upper_bound=0)
            # solve the model
            sol_pam = pamodel.optimize()
            if pamodel.solver.status == "optimal" and pamodel.objective.value != 0:
                # substrate_range += [abs(substrate)]
                substrate_range[abs(substrate)] = [abs(pamodel.reactions.get_by_id(substrate_uptake_id).flux)]
                fluxes.append(sol_pam.fluxes)
                if save_fluxes_esc and sensitivity: self.save_pamodel_simulation_results(substrate_uptake_rate=substrate,
                                                                         substrate_uptake_reaction=substrate_uptake_id,
                                                                         bin_id= "no bins")

            # reset substrate_uptake_rate
            pamodel.change_reaction_bounds(substrate_uptake_id,
                                              lower_bound=previous_bounds[0], upper_bound=previous_bounds[1])

        return fluxes, substrate_range
