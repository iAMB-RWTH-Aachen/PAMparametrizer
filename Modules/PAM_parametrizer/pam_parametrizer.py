from typing import Union, Tuple, Iterable
import numpy as np
import os
import time
import matplotlib.pyplot as plt
import matplotlib as mpltlib
from matplotlib.colors import to_hex
from PAModelpy.PAModel import PAModel
import pandas as pd
from scipy.stats import linregress
from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import StandardScaler
import re
from cobra import DictList

from .PAM_data_classes import ValidationData, HyperParameters, ParametrizationResults
from Modules.genetic_algorithm_parametrization import GAPOGaussian as GAPOGauss
from Modules.utils.error_calculation import calculate_r_squared_for_reaction, nanaverage
from Modules.utils.sampling_functions import adaptive_sampling
from Modules.utils.sector_config_functions import get_model_simulations_vs_sector, perform_linear_regression, change_translational_sector_with_config_dict


class PAMParametrizer():

    #from Schmidt et al(2016):
    TRANSL_SECTOR_INTERCEPT_vs_MU = 0.04738115630907698#g / g_cdw
    TRANSL_SECTOR_SLOPE_vs_MU = 0.04806975534478209 #g / g_cdw / h
    TRANSL_SECTOR_INTERCEPT_vs_GLC = 0.046136644909661115#g / g_cdw
    TRANSL_SECTOR_SLOPE_vs_GLC = -0.004340256958025938 #g / mmol_glc / h
    MEASURED_PROTEIN_FRACTION = 0.55*0.55

    def __init__(self, pamodel:PAModel,
                 validation_data: Union[DictList[ValidationData], list, ValidationData],
                 hyperparameters: HyperParameters = HyperParameters(),
                 substrate_uptake_id: str = "EX_glc__D_e",
                 max_substrate_uptake_rate: Union[float, int] = 0,
                 min_substrate_uptake_rate: Union[float, int] = -11,
                 sensitivity: bool = True,
                 enzymes_to_evaluate:list = []):

        self.core_genetic_algorithm = None
        self.pamodel = pamodel
        #change total protein constraint to equality constraint for better fittinh
        self._set_total_protein_constraint_to_equality()
        if not hasattr(validation_data, "__iter__"): validation_data = [validation_data]
        self.validation_data = DictList(validation_data)
        self.hyperparameters = hyperparameters

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
        if not sensitivity:
            self.enzymes_to_evaluate = enzymes_to_evaluate

        # attributes for keeping track of the workflow
        self.iteration = 0
        self.bins = list(range(self.hyperparameters.number_of_bins))
        self.result_figure_file = os.path.join(
            "Results", f"pam_parametrizer_progress_{hyperparameters.filename_extension}.png")
        self.result_diagnostics_file = os.path.join(
            "Results", f"pam_parametrizer_diagnostics_{hyperparameters.filename_extension}.xlsx")


    def run(self, remove_subruns:bool = True, binned:str = "all") -> None:
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
        self.calculate_translational_sector_for_multiple_csources()

        #setup plot to visualize progress
        fig, axs = self.plot_valid_data()
        fig = self.plot_simulation(fig=fig, axs=axs, color="#010328")

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

            if self.evaluate_and_save_results_of_iteration(start_time_iteration, files_to_remove, remove_subruns, fig, axs) is not None:
                #if the error is converging, optimize random enzyme parameters
                self.iteration +=1
                files_to_remove = self.perform_iteration_without_bins(random=True)
                self.evaluate_and_save_results_of_iteration(start, files_to_remove, remove_subruns, fig, axs)


        self.save_final_diagnostics(figure = fig)
        plt.close(fig)

    def calculate_translational_sector_for_multiple_csources(self):
        #generate a pam with only the translational sector
        pamtransl = self.pamodel.copy()
        pamtransl.sensitivity = False #this will speed the simulations up
        #remove total protein to remove protein relations
        pamtransl.remove_cons_vars([pamtransl.constraints[pamtransl.TOTAL_PROTEIN_CONSTRAINT_ID]])

        for vd in self.validation_data:
            if vd.translational_sector_config is not None: continue
            sub_upt_id = vd.id
            #to prevent overflow metabolism, only select the lower growth rates to derive the equation
            substrate_range = self._get_substrate_range_lower_substrate_conc(vd.validation_range)
            #get the simulations for relatively low substrate uptake rates
            simulation_results = get_model_simulations_vs_sector(pamodel = pamtransl,
                                                                 sub_uptake_rxn = sub_upt_id,
                                                                 rxn_id_to_relate_to = self.pamodel.BIOMASS_REACTION,
                                                                 substrate_range = substrate_range,
                                                                 intercept = self.TRANSL_SECTOR_INTERCEPT_vs_MU,
                                                                 slope = self.TRANSL_SECTOR_SLOPE_vs_MU)
            if len(simulation_results)==0 or len(simulation_results[sub_upt_id].unique()) ==1: # model was infeasible for all datapoints or did not change at all
                slope, intercept = self.TRANSL_SECTOR_SLOPE_vs_GLC, self.TRANSL_SECTOR_INTERCEPT_vs_GLC
            else:
                slope, intercept = perform_linear_regression(
                x=simulation_results[sub_upt_id], y=simulation_results["translational_protein"])
                # slope, intercept = slope * self.MEASURED_PROTEIN_FRACTION, intercept*self.MEASURED_PROTEIN_FRACTION
            vd.translational_sector_config = {"slope":slope,"intercept":intercept}


    def perform_iteration_in_bins(self, start_time:float) -> list:
        self._init_results_objects()
        # 1. Get the binned substrate values and adjust binsize if required
        binned_substrates = self._bin_substrate_uptake_rates()


        # 2. Run model in bins, get sensitivities and calculate errors
        for subst_uptake_id, binned_substrate in binned_substrates.items():
            self._change_translational_sector_for_substrate(substrate_uptake_id=subst_uptake_id)
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
            esc_df = self.pamodel.enzyme_sensitivity_coefficients
            esc_df["rxn_id"] = esc_df["rxn_id"].str.split(",")
            esc_df = esc_df.explode("rxn_id", ignore_index=True)
            enzymes_to_evaluate = self._parse_enzymes_to_evaluate(
                esc_df.sample(
                    self.hyperparameters.number_of_kcats_to_mutate
                ).rename({"coefficient": "mean"}, axis=1)
            )

        # 3. run genetic algorithm for a range of substrate uptake rates
        self.run_genetic_algorithm(enzymes_to_evaluate=enzymes_to_evaluate,
                                   filename_extension=f"final_run_{self.iteration}")

        # 3. Restart genetic algorithm using populations from bins
        files_to_remove = self._get_genetic_algorithm_json_files()


        return files_to_remove

    def evaluate_and_save_results_of_iteration(self, start_time_iteration:float, files_to_remove:list,
                                               remove_subruns:bool, fig: plt.Figure, axs: plt.Axes):
        self.reparametrize_pam()
        if self._pamodel_is_feasible:
            self._init_results_objects()
            # visualize results
            fig = self.plot_simulation(fig, axs, save_esc=True)

            for valid_data in self.validation_data:
                sampled_data = valid_data.sampled_valid_data
                substrate_uptake_id = valid_data.id
                substrate_rates_for_error = sampled_data[substrate_uptake_id+"_ub"].to_list()

                fluxes, substrate_rates = self.run_simulations_to_plot(substrate_uptake_id,
                                                                       substrate_rates=substrate_rates_for_error)

                for simulation_result, substrate_rate in zip(fluxes,substrate_rates):
                    self.parametrization_results.add_fluxes_from_fluxdict(flux_dict=simulation_result,
                                                                          bin_id="final",
                                                                          substrate_reaction_id= substrate_uptake_id,
                                                                          substrate_uptake_rate=substrate_rate,
                                                                          fluxes_abs=False)
            self.final_error = self.calculate_final_error()


            self.save_diagnostics(computational_time=time.perf_counter() - start_time_iteration)
            # 6. Display progress and repeat
            print("time elapsed: ", time.perf_counter() - start_time_iteration, "sec, ",
                  (time.perf_counter() - start_time_iteration) / 60, "min",
                  (time.perf_counter() - start_time_iteration) / 3600, "hour")
            print("Done with iteration number ", self.iteration)
            print("-------------------------------------------------------------------------------------------")

            if remove_subruns: self._remove_result_files(files_to_remove)
            if self._error_is_converging():
                print("Stopped simulations because error is converging")
                return np.nan
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
                                                bin_id = bin_id, substrate_uptake_rates =substrate_rates)
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
                                   filename_extension=f"iteration_{self.iteration}_bin_{bin_id}")

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
        self._change_translational_sector_for_substrate(substrate_uptake_reaction)
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

        print("Substrate uptake rate ", substrate_uptake_rate, " mmol/gcdw/h")
        with self.pamodel:
            # change glucose uptake rate
            if substrate_uptake_rate < 0:
                self.pamodel.change_reaction_bounds(rxn_id=substrate_uptake_reaction,
                                                    lower_bound=substrate_uptake_rate, upper_bound=0)
            else:
                self.pamodel.change_reaction_bounds(rxn_id=substrate_uptake_reaction,
                                                lower_bound=0, upper_bound=substrate_uptake_rate)
            # solve the model
            self.pamodel.optimize()
            if self.pamodel.solver.status == "optimal" and self.pamodel.objective.value != 0:
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

        self.parametrization_results.add_fluxes(self.pamodel, bin_id, substrate_uptake_reaction, substrate_uptake_rate,
                                                fluxes_abs=False)
        self.parametrization_results.add_enzyme_sensitivity_coefficients(self.pamodel.enzyme_sensitivity_coefficients,
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
        validation_df = valid_data.sampled_valid_data
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
            substrate_uptake_id = valid_data.id
            reactions_to_validate = valid_data._reactions_to_validate
            validation_df = valid_data.sampled_valid_data#.apply(lambda x: x.abs() if x.dtype!='object' else x)

            # self._init_validation_df(bin_information=[self.min_substrate_uptake_rate, self.max_substrate_uptake_rate])
            error += [nanaverage(self._calculate_error_for_reactions(substrate_uptake_id = substrate_uptake_id,
                                                         validation_df=validation_df,
                                                         reactions_to_validate= reactions_to_validate,
                                                         bin_id="final"))]


            #remove the simulations from the result dataframe
            self.parametrization_results.remove_simulations_from_flux_df(substrate_uptake_id,"final")
        final_error = np.nanmean(error)
        print("The final error is ", final_error)

        return final_error

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

    def run_genetic_algorithm(self, enzymes_to_evaluate:dict, filename_extension:str) -> None:
        """ Function to run the genetic algorithm.

        Results will be stored in a json, pickle and xlsx file with
        filename defined by: self.hyperparameters.genetic_algorithm_filename_base + filename_extension

        Args:
            enzymes_to_evaluate (dict): Dictionary with required information about which enzyme parameters should be optimized.
                Format:
                {
                    'enzyme_id': {
                        'reaction': 'reaction_id',
                        'kcats': {'f':kcat_value_forward, 'b':kcat_value_reverse},
                        'sensitivity': esc_value
                    }
                }
            filename_extension: extension of base filename to save the result of the genetic algorithm
        """
        substrate_rates = dict()
        translation_configs = dict()
        for valid_data in self.validation_data:
            #need to make sure that there is data sampled in the range to evaluate for this specific substrate
            if valid_data.sampled_valid_data is not None:
                substrate_rates[valid_data.id] = valid_data.sampled_valid_data[valid_data.id+"_ub"]
            translation_configs[valid_data.id] = valid_data.translational_sector_config

        ga = self._init_genetic_algorithm(substrate_uptake_rates=substrate_rates,
                                            enzymes_to_evaluate=enzymes_to_evaluate,
                                            translational_sector_config=translation_configs,
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
        translational_sector_config = {vd.id: vd.translational_sector_config for vd in self.validation_data}

        ga = self._init_genetic_algorithm(substrate_uptake_rates= substrate_rates,
                                          enzymes_to_evaluate= enzymes_to_evaluate,
                                          translational_sector_config= translational_sector_config,
                                          filename_extension= f"final_run_{self.iteration}")
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

        print("after reparametrization the model is feasible", self._pamodel_is_feasible())


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

    def save_final_diagnostics(self, figure:plt.Figure) -> None:
        """Save the results of all iterations to excel and png

        Args:
            figure: figure to save (progress figure)
        """
        figure.savefig(self.result_figure_file, dpi=100, bbox_inches="tight")
        transl_sector_config = pd.DataFrame(columns = ["substrate_uptake_id", "slope", "intercept"])
        for vd in self.validation_data:
            transl_sector_config.loc[len(transl_sector_config)] = [vd.id, vd.translational_sector_config["slope"],
                                                                   vd.translational_sector_config["intercept"]]

        with pd.ExcelWriter(self.result_diagnostics_file) as writer:
            # Write each DataFrame to a specific sheet
            self.parametrization_results.best_individuals.to_excel(writer, sheet_name="Best_Individuals", index=False)
            self.parametrization_results.computational_time.to_excel(writer, sheet_name="Computational_Time",
                                                                     index=False)
            self.parametrization_results.final_errors.to_excel(writer, sheet_name="Final_Errors", index=False)
            transl_sector_config.to_excel(writer, sheet_name="translational_sector", index = False)

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


    ###########################################################################################################
    #WORKER FUNCTIONS
    ###########################################################################################################
    def _change_translational_sector_for_substrate(self, substrate_uptake_id: str) -> None:
        transl_sector_config = self.validation_data.get_by_id(substrate_uptake_id).translational_sector_config
        if transl_sector_config is None: return # use default if other parameterization is not provided
        change_translational_sector_with_config_dict(self.pamodel, transl_sector_config, substrate_uptake_id)

    def _get_substrate_range_lower_substrate_conc(self, validation_range:list[Union[int, float]], number_of_steps: int = 5) -> Iterable:
        #only loop over the low growth rates, to prevent overflow like metabolism to interfere with the derivation of a linear equation
        # Determine the range to loop over based on absolute values
        start = min(abs(validation_range[0]), abs(validation_range[1])) / 2
        end = max(abs(validation_range[0]), abs(validation_range[1])) / 2
        step = abs(start-end)/number_of_steps

        # Create the substrate range
        substrate_range = np.arange(start, end, step)

        # Adjust the sign if the original range was negative
        if validation_range[0] < 0 and validation_range[1] < 0:
            substrate_range = -substrate_range[::-1]
        return substrate_range

    def _pamodel_is_feasible(self):
        #check for feasibility at a mid substrate uptake rate
        self.pamodel.test((self.max_substrate_uptake_rate -self.min_substrate_uptake_rate/2))
        return self.pamodel.solver.status == "optimal"

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

    def _init_genetic_algorithm(self, substrate_uptake_rates: dict,
                                enzymes_to_evaluate: dict,
                                translational_sector_config: dict,
                                filename_extension:str) -> GAPOGauss:
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
                    'enzyme_id': {
                        'reaction': 'reaction_id',
                        'kcats': {'f':kcat_value_forward, 'b':kcat_value_reverse},
                        'sensitivity': esc_value
                    }
                }
            translational_sector_config (dict of dict): Dictionary with the slope and intercept of the translational
                sector configuration for each substrate.
                Format:
                {'substrate_uptake_id':{
                    'slope':float, #slope in g/mmol/h
                    'intercept':float #intercept in g/mmol
                    }
                }
            filename_extension (str): Extension of the base filename to save the result of the genetic algorithm.

        Returns:
            ga: Core genetic algorithm object ready to be started.
        """

        results_filename = self.hyperparameters.genetic_algorithm_filename_base + filename_extension

        ga = self.hyperparameters.genetic_algorithm(
            model=self.pamodel,
            enzymes_to_eval= enzymes_to_evaluate,
            translational_sector_config = translational_sector_config,
            valid_data = {valid_data.id: valid_data.sampled_valid_data for valid_data in self.validation_data if valid_data.sampled_valid_data is not None},
            filename_save = results_filename,
            substrate_uptake_id = self.substrate_uptake_id,
            substrate_uptake_rates = substrate_uptake_rates, # bin_info: [start, stop, step]
            objective_id = self.validation_data[0]._get_biomass_reactions(),
            **self.hyperparameters.genetic_algorithm_hyperparams
        )
        return ga

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

    def _calculate_error_for_reactions(self, substrate_uptake_id:str,
                                       validation_df: pd.DataFrame,
                                       reactions_to_validate: list,
                                        bin_id: Union[float, int] = None) -> float:
        # calculate error for different exchange rates
        error = []
        flux_df = self.parametrization_results.flux_results.get_by_id(substrate_uptake_id).fluxes_df

        if len(flux_df) == 0:  # means model is infeasible
            return -1
        # check if we want to calculate the error for a single bin
        if bin_id is not None: flux_df = flux_df[flux_df["bin"] == bin_id]


        # if bin_id is None: bin_id = "no bins"
        for rxn in reactions_to_validate: #+ self.validation_data._get_biomass_reactions():
            # only select the rows which are filled with data
            validation_data = validation_df.dropna(axis=0, subset=rxn)
            # if there are no reference data points, continue to the next reaction
            if len(validation_data) == 0:
                error += [np.nan]
                continue

            r_squared = calculate_r_squared_for_reaction(rxn, validation_df, substrate_uptake_id,
                                                               flux_df)

            error += [r_squared]
        return error


    def _parse_enzymes_to_evaluate(self, esc_topn_df: pd.DataFrame) -> dict:
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

        enzymes_to_evaluate = {}
        for index, row in esc_topn_df.iterrows():
            enzyme_id = row["enzyme_id"]
            rxn_id = row["rxn_id"]
            enzyme_dict = {}
            enzyme_dict["sensitivity"] = row["mean"]
            #create the connection to the catalytic reaction related to the enzyme
            if "CE_" not in rxn_id: rxn_id = f"CE_{rxn_id}_{enzyme_id}"
            enzyme_dict["reaction"] = rxn_id
            kcat_dict = self.pamodel.enzymes.get_by_id(enzyme_id).rxn2kcat[rxn_id]
            enzyme_dict["kcats"] = {dir: 1/(kcat* 3600 * 1e-6) for dir, kcat in kcat_dict.items() if kcat>0}
            enzymes_to_evaluate[enzyme_id] = enzyme_dict
        return enzymes_to_evaluate

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
        if len(esc_results_df) == 0: return None
        esc_topn_df = self._select_topn_enzymes(esc_results_df,
                                                nmbr_kcats_to_pick)
        enzymes_to_evaluate = self._parse_enzymes_to_evaluate(esc_topn_df)

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
        self.pamodel.change_kcat_value(enzyme_id=enzyme_id, kcats=kcat_dict)
        for rxn, kcat in kcat_dict.items():
            for direction, kcat in kcat.items():
                self.pamodel._change_kcat_in_enzyme_constraint(rxn, enzyme_id, direction, kcat)


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
        self.pamodel.constraints[self.pamodel.TOTAL_PROTEIN_CONSTRAINT_ID].lb = self.pamodel.constraints[self.pamodel.TOTAL_PROTEIN_CONSTRAINT_ID].ub

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
            x = [abs(glc) for glc in valid_data.valid_data[self.substrate_uptake_id +"_ub"]]
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

        #TODO: set bounds dynamically
        # axs.flatten()[-1].plot([-10, 10], [-10,10], linestyle = "dotted")
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
                        color:int= None) -> plt.Figure:

        if color is None:
            #adjust color to visualize progress
            #get viridis color palette
            cmap = plt.get_cmap("viridis")
            color = to_hex(cmap(self.iteration/(self.hyperparameters.threshold_iteration)))
            norm = mpltlib.colors.BoundaryNorm(list(range(self.hyperparameters.threshold_iteration+1)), cmap.N)

        fluxes_dict = {}
        substrate_range_dict = {}

        for substrate_id in self.substrate_uptake_ids:
            if substrate_id == self.substrate_uptake_id: substrate_range = None
            else: substrate_range = self.validation_data.get_by_id(substrate_id).valid_data[substrate_id]

            fluxes, substrate_range = self.run_simulations_to_plot(substrate_uptake_id=substrate_id,
                                                                   substrate_rates=substrate_range,
                                                                    save_fluxes_esc=save_esc)
            if len(fluxes) > 0: # only plot feasible model results
                fluxes_dict[substrate_id] = fluxes
                substrate_range_dict[substrate_id] = substrate_range


        alpha = 1
        #only plot detailed info for the 'main' substrate
        for r, ax in zip(self.validation_data.get_by_id(self.substrate_uptake_id)._reactions_to_plot, axs.flatten()[:-1]):
            # plot data
            line = ax.plot(substrate_range_dict[self.substrate_uptake_id], [abs(f[r]) for f in fluxes_dict[self.substrate_uptake_id]], linewidth=2.5,
                           zorder=5, color=color, alpha=alpha)
        # Plot a flux comparison to experimental data for the other fluxes
        for substr_id, fluxes in fluxes_dict.items():
            valid_data = self.validation_data.get_by_id(substr_id)
            if substr_id != self.substrate_uptake_id and valid_data.sampled_valid_data is not None:
                for reaction in valid_data._reactions_to_validate:
                    exp_measurements = valid_data.sampled_valid_data[reaction]
                    simulations = [f[reaction] for f in fluxes]
                    axs.flatten()[-1].scatter(exp_measurements, simulations, color = color, alpha = alpha)

        # Add colorbar
        if self.iteration == 1:
            cbar = fig.colorbar(mpltlib.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axs)
            cbar.set_label("Iteration")

        fig.canvas.draw()
        fig.canvas.flush_events()
        if return_fluxes: return fig, fluxes_dict, substrate_range_dict
        return fig

    def run_simulations_to_plot(self, substrate_uptake_id:str,
                                substrate_rates: Union[np.array, list, pd.Series] = None,
                                save_fluxes_esc:bool = False) -> Tuple[list, list]:
        fluxes = list()
        substrate_range = list()
        #get the old bounds for resetting

        lb, ub = -self.pamodel.constraints[substrate_uptake_id+ "_lb"].ub, self.pamodel.constraints[substrate_uptake_id+ "_ub"].ub
        self._change_translational_sector_for_substrate(substrate_uptake_id)

        if substrate_rates is None:
            step = (self.max_substrate_uptake_rate-self.min_substrate_uptake_rate)/10
            substrate_rates = np.arange(self.min_substrate_uptake_rate, self.max_substrate_uptake_rate, step)
        for substrate in substrate_rates:

            if substrate>=0:
                self.pamodel.change_reaction_bounds(rxn_id=substrate_uptake_id,
                                                lower_bound=0, upper_bound=substrate)
            else:
                self.pamodel.change_reaction_bounds(rxn_id=substrate_uptake_id,
                                                    lower_bound=substrate, upper_bound=0)
            # solve the model
            sol_pam = self.pamodel.optimize()
            if self.pamodel.solver.status == "optimal" and self.pamodel.objective.value != 0:
                substrate_range += [abs(substrate)]
                fluxes.append(sol_pam.fluxes)
                if save_fluxes_esc: self.save_pamodel_simulation_results(substrate_uptake_rate=substrate,
                                                                         substrate_uptake_reaction=substrate_uptake_id,
                                                                         bin_id= "no bins")


            # reset substrate_uptake_rate
            if substrate<0:
                self.pamodel.change_reaction_bounds(substrate_uptake_id,
                                              lower_bound=0, upper_bound=1e6)
            else:
                self.pamodel.change_reaction_bounds(substrate_uptake_id,
                                              lower_bound=-1e6, upper_bound=0)

            # if substrate >= 0:
            #     self.pamodel.change_reaction_bounds(rxn_id=substrate_uptake_id,
            #                                         lower_bound=lb, upper_bound=0)
            # else:
            #     self.pamodel.change_reaction_bounds(rxn_id=substrate_uptake_id,
            #                                         lower_bound=0, upper_bound=ub)

        return fluxes, substrate_range

    def cluster_parametrization_results(self, parameter_sets_list, fitness_values_list, n_clusters):
        """
        Clusters the concatenated parameter sets obtained from multiple parametrization runs,
        considering both the parameter values and the enzyme IDs, and estimates the fitness value for each cluster.

        Args:
        - parameter_sets_list (list of lists of tuples): List of parameter sets obtained from each parametrization run,
                                                         where each parameter set is represented as a tuple
                                                         (enzyme_id, parameter_set).
        - fitness_values_list (list of lists): List of fitness values obtained from each parametrization run.
        - n_clusters (int): Number of clusters to form.

        Returns:
        - clustered_results (list of tuples): Clustered parameter sets, where each tuple contains
                                               the cluster label, parameter sets, enzyme IDs, and estimated fitness value.
        """
        # Concatenate parameter sets and enzyme IDs from all runs
        concatenated_parameter_sets = []
        for run in parameter_sets_list:
            concatenated_parameter_sets.extend(run)

        # Extract parameter sets and enzyme IDs
        parameter_sets = np.array([param_set for enzyme_id, param_set in concatenated_parameter_sets])
        enzyme_ids = np.array([enzyme_id for enzyme_id, param_set in concatenated_parameter_sets])

        # Standardize parameter sets (optional but recommended for clustering)
        scaler = StandardScaler()
        parameter_sets_standardized = scaler.fit_transform(parameter_sets)

        # Concatenate parameter sets and enzyme IDs horizontally
        data = np.hstack((parameter_sets_standardized, enzyme_ids.reshape(-1, 1)))

        # Initialize Agglomerative Clustering with Ward linkage
        clustering = AgglomerativeClustering(n_clusters=n_clusters, linkage="ward")

        # Fit clustering model
        cluster_labels = clustering.fit_predict(data)

        # Compute mean fitness value for each cluster
        clustered_fitness_values = []
        for i in range(n_clusters):
            cluster_indices = np.where(cluster_labels == i)[0]
            fitness_values = [fitness_values_list[idx] for idx in cluster_indices]
            mean_fitness = np.mean(fitness_values)
            clustered_fitness_values.append(mean_fitness)

        # Clustered results (cluster label, parameter sets, enzyme IDs, and estimated fitness value)
        clustered_results = []
        for label, cluster_fitness in enumerate(clustered_fitness_values):
            cluster_indices = np.where(cluster_labels == label)[0]
            cluster_parameter_sets = [parameter_sets[i] for i in cluster_indices]
            cluster_enzyme_ids = [enzyme_ids[i] for i in cluster_indices]
            clustered_results.append((label, cluster_parameter_sets, cluster_enzyme_ids, cluster_fitness))

        return clustered_results