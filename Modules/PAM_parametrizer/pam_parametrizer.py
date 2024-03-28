from typing import Union
import numpy as np
import os
import time
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from PAModelpy.PAModel import PAModel
import pandas as pd
from scipy.stats import linregress
from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import StandardScaler
import re

from .PAM_data_classes import ValidationData, HyperParameters, ParametrizationResults
from Modules.genetic_algorithm_parametrization import GAPOGaussian as GAPOGauss


class PAMParametrizer():
    def __init__(self, pamodel:PAModel,
                 validation_data: ValidationData = ValidationData,
                 hyperparameters: HyperParameters = HyperParameters,
                 substrate_uptake_id: str = 'EX_glc__D_e',
                 max_substrate_uptake_rate: Union[float, int] = 11,
                 min_substrate_uptake_rate: Union[float, int] = 0,
                 sensitivity: bool = True,
                 enzymes_to_evaluate:list = []):

        self.core_genetic_algorithm = None
        self.pamodel = pamodel
        #change total protein constraint to equality constraint for better fittinh
        self._set_total_protein_constraint_to_equality()
        self.validation_data = validation_data
        self.hyperparameters = hyperparameters
        self.parametrization_results = ParametrizationResults()
        self.enzyme_ids = [enzyme.id for enzyme in pamodel.enzymes]

        self.substrate_uptake_id = substrate_uptake_id
        self.min_substrate_uptake_rate = min_substrate_uptake_rate
        self.max_substrate_uptake_rate = max_substrate_uptake_rate
        self.fluxes_df = pd.DataFrame()
        self.final_error = 0
        self.error_fraction_to_deviate_between_runs = 0.1

        self.sensitivity = sensitivity
        if not sensitivity:
            self.enzymes_to_evaluate = enzymes_to_evaluate

        # attributes for keeping track of the workflow
        self.iteration = 0
        self.bins = list(range(self.hyperparameters.number_of_bins))
        self.result_figure_file = os.path.join(
            'Results', f'pam_parametrizer_progress_{hyperparameters.filename_extension}.png')
        self.result_diagnostics_file = os.path.join(
            'Results', f'pam_parametrizer_diagnostics_{hyperparameters.filename_extension}.xlsx')


    def run(self, remove_subruns:bool = True, binned:str = 'all') -> None:
        """ Run the parametrization framework.

        For each iteration of parametrization the following steps are taken:
        0. Initiate objects which are needed to store the results
        1. Split the range of substrate uptake rates to evaluate in a number of bins
        2. For each bin:
            i. Run simulations and calculate error to experimental measurements
            ii. Select the parameters to optimize
            iii. Run a genetic algorithm to optimize the selected parameters
        3. Select the enzymes to evaluate for the entire range of substrate uptake rates
        4. Use the results of the previous genetic algorithm runs to initialize a new genetic algorithm
        5. Run the genetic algorithm and use the best parameter set to reparametrize the model

        This is repeated until a max number of iterations is reached or until the simulation error is below a threshold

        Args:
            remove_subruns: If True, the results from the genetic algorithm runs within the bins are removed.
            The results of the final run won't be removed.
        """
        #setup plot to visualize progress
        fig, axs = self.plot_valid_data()
        fig = self.plot_simulation(fig=fig, axs=axs, color='#010328')

        #keep track of time for computational performance
        start = time.perf_counter()
        if binned == 'before':
            self.validation_data.sampled_valid_data_df = self._adaptive_sampling(self.validation_data.valid_data_df)
            self._init_results_objects()
            files_to_remove = self.perform_iteration_in_bins(start)
            self.evaluate_and_save_results_of_iteration(start, files_to_remove, remove_subruns, fig, axs)

        #perform parametrization until a max number of iterations
        while (self.iteration+1 <= self.hyperparameters.threshold_iteration) & (self.final_error <= self.hyperparameters.threshold_error):
            # 0. initiate structures and objects to save results and keep track of process
            start_time_iteration = time.perf_counter()
            self.iteration += 1

            if binned == 'all': #binning in all iterations
                # 1-4. Run simulations and genetic algorithm in bins and
                # restart a genetic algorithm using the resulting populations
                files_to_remove = self.perform_iteration_in_bins(start)
            else: # no binning
                # 1-4. Run genetic algorithm without bins
                files_to_remove = self.perform_iteration_without_bins()

            # 5. Reparametrize the model with best results and save/visualize results
            self.evaluate_and_save_results_of_iteration(start_time_iteration, files_to_remove, remove_subruns, fig, axs)

            # 6. Display progress and repeat
            print('time elapsed: ', time.perf_counter() - start, 'sec, ', (time.perf_counter() - start) / 60, 'min',
                  (time.perf_counter() - start) / 3600, 'hour')
            print('Done with iteration number ', self.iteration)
            print('-------------------------------------------------------------------------------------------')

        self.save_final_diagnostics(figure = fig)
        plt.close(fig)

    def perform_iteration_in_bins(self, start_time:float) -> list:
        self._init_results_objects()
        # 1. Get the binned substrate values and adjust binsize if required
        binned_substrate = self._bin_substrate_uptake_rates()

        # 2. Run model in bins, get sensitivities and calculate errors
        for bin_id, bin_info in binned_substrate.items():
            self.process_bin(bin_info, bin_id)
            # print running time to check on progress
            print('time elapsed: ', time.perf_counter() - start_time, 'sec, ', (time.perf_counter() - start_time) / 60, 'min',
                  (time.perf_counter() - start_time) / 3600, 'hour\n')


        # 3. Restart genetic algorithm using populations from bins
        files_to_remove = self.restart_genetic_algorithm()

        return files_to_remove

    def perform_iteration_without_bins(self) -> list:
        # 1. If there are no results yet, perform simulations for a range of glc_uptake_rates
        if self.iteration == 1:
            self._init_results_objects()
            self.run_simulations_to_plot(save_fluxes_esc=True)

        # 2. Get the most sensitive enzymes from all runs
        if self.sensitivity:
            enzymes_to_evaluate = self._determine_enzymes_to_evaluate_for_all_bins(
                nmbr_kcats_to_pick=self.hyperparameters.number_of_kcats_to_mutate)
        else:
            enzymes_to_evaluate = self.enzymes_to_evaluate

        # 3. run genetic algorithm for a range of substrate uptake rates
        step = (self.max_substrate_uptake_rate - self.min_substrate_uptake_rate) / 10
        substrate_rates = list(np.arange(self.min_substrate_uptake_rate, self.max_substrate_uptake_rate, step))
        ga = self._init_genetic_algorithm(substrate_uptake_rates=substrate_rates,
                                          enzymes_to_evaluate=enzymes_to_evaluate,
                                          filename_extension=f'final_run_{self.iteration}')
        ga.start()

        # 3. Restart genetic algorithm using populations from bins
        files_to_remove = self._get_genetic_algorithm_json_files()

        return files_to_remove

    def evaluate_and_save_results_of_iteration(self, start_time_iteration:float, files_to_remove:list,
                                               remove_subruns:bool, fig: plt.Figure, axs: plt.Axes):
        self.reparametrize_pam()
        if self._pamodel_is_feasible:
            self._init_results_objects()
            # visualize results
            fig, fluxes, substrate_rates = self.plot_simulation(fig, axs, return_fluxes=True, save_esc=True)
            self.final_error = self.calculate_final_error(fluxes, substrate_rates)
            self.save_diagnostics(computational_time=time.perf_counter() - start_time_iteration)
            if remove_subruns: self._remove_result_files(files_to_remove)
        else:
            print('Newly parametrized model is not feasible')
            return

    def process_bin(self, bin_information:list, bin_id: str) -> None:
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
            f'The following range of substrate uptake rates will be analyzed: {bin_information[0]} - {bin_information[1]}')
        if self._init_validation_df(bin_information) is None: return
        self.run_pamodel_simulations_in_bin(bin_information={bin_id: bin_information})
        # calculate the error for the different exchange rates
        self.calculate_error(bin_information, bin_id)
        # if there is no solution, continue to the next bin
        if len(self.parametrization_results.esc_df[self.parametrization_results.esc_df['bin'] == bin_id]) == 0:
            return
        # get the most sensitive enzymes
        top_enzyme_sensitivities = self.determine_most_sensitive_enzymes(bin_id,
                                                                         self.hyperparameters.number_of_kcats_to_mutate)
        self.determine_bin_to_split(top_enzyme_sensitivities, bin_id)

        self.run_genetic_algorithm(bin_info=bin_information,
                                   esc_topn_df=top_enzyme_sensitivities,
                                   filename_extension=f'iteration_{self.iteration}_bin_{bin_id}')

    def run_pamodel_simulations_in_bin(self, bin_information:dict) -> None:
        """ Run pamodel in range of substrate uptake rates defined in bin.

        Use the range of substrate uptake rate indicated in the bin_information to run simulations with the PAModel.
        Saves the simulation results in the parametrization_results data class.

        Args:
            bin_information: dictionary with bin_id:[start, stop, step] key:value pairs, where start, stop and step
                            relate to the start, end and stepsize of the substrate uptake rate range
        """
        for bin_id, bin_info in bin_information.items():
            start, stop, step = bin_info[0], bin_info[1], bin_info[2]
            substrate_range = self.validation_data.valid_data_df[self.substrate_uptake_id+'_ub']
            print(
                f'The following range of substrate uptake rates will be analyzed: {start - stop}')
            for substrate_uptake_rate in substrate_range:
                self.run_pamodel_simulation(substrate_uptake_rate, bin_id)

    def run_pamodel_simulation(self, substrate_uptake_rate: Union[float, int],
                               bin_id: Union[str, float, int]) -> None:
        """ Running pamodel simulation for specific substrate uptake rate.

        Running PAModel simulations and saving the resulting fluxes and enzymes sensitivities coefficient
        for later analysis

        Args:
            substrate_uptake_rate: used to constrain the PAModel
            bin_id: identifier of the bin in which this simulation is run (for saving purposes)
        """

        print('Substrate uptake rate ', substrate_uptake_rate, ' mmol/gcdw/h')
        with self.pamodel:
            # change glucose uptake rate
            if substrate_uptake_rate < 0:
                self.pamodel.change_reaction_bounds(rxn_id=self.substrate_uptake_id,
                                                    lower_bound=substrate_uptake_rate, upper_bound=0)
            else:
                self.pamodel.change_reaction_bounds(rxn_id=self.substrate_uptake_id,
                                                lower_bound=0, upper_bound=substrate_uptake_rate)
            # solve the model
            self.pamodel.optimize()
            if self.pamodel.solver.status == 'optimal' and self.pamodel.objective.value != 0:
                self.save_pamodel_simulation_results(abs(substrate_uptake_rate), bin_id)

    def save_pamodel_simulation_results(self, substrate_uptake_rate: Union[float, int],
                                        bin_id: Union[str, float, int]) -> None:
        """ Saving simulation results.

            Saving the resulting fluxes and enzymes sensitivities coefficient from a successful PAModel simulation
            for later analysis

            Args:
                substrate_uptake_rate: used to constrain the PAModel
                bin_id: identifier of the bin in which this simulation is run (for saving purposes)
        """

        self.parametrization_results.substrate_range += [substrate_uptake_rate]
        self.parametrization_results.add_fluxes(self.pamodel, bin_id, substrate_uptake_rate)
        self.parametrization_results.add_enzyme_sensitivity_coefficients(self.pamodel.enzyme_sensitivity_coefficients,
                                                                          bin_id, substrate_uptake_rate)

    def calculate_error(self, bin_information: list, bin_id: Union[float, int, str]) -> None:
        """ Evaluate the model simulations compared to a reference dataset within a specified substrate uptake range.

        This method calculates the average difference between model simulations and reference data for the total substrate
        uptake range and available reactions within the specified bin.

        Notes:
        - The reference dataset is filtered to include only data points within the specified substrate uptake range.
        - Errors are calculated for different exchange rates and biomass reactions.
        - The '_calculate_r_squared_for_reaction' method is used to determine the R-squared value for each reaction.
        - The calculated errors are stored in the 'error_df' attribute for further analysis.

        Args:
            bin_information: List containing upper and lower bounds of the substrate uptake range for the bin.
            bin_id: Identifier for the bin under consideration.
        """
        # calculate error for different exchange rates
        error = self._calculate_error_for_reactions(validation_df = self.validation_data.sampled_valid_data_df,
                                                    # substrate_range = self.validation_data.sampled_valid_data_df[self.substrate_uptake_id+'_ub'],
                                                    bin_id = bin_id)

        error_df = self.parametrization_results.error_df
        self.parametrization_results.error_df.loc[len(error_df)] = [bin_id] + error

    def calculate_final_error(self, fluxes: list, substrate_rates:list) -> float:
        """ Calculate the simulation error of the reparametrized model

        Args:
            fluxes: a list of flux solutions (pd.Series) as returned by the cobra model.optimize() function.
            substrate_rates: the list of substrate uptake rates which were used to constrain the simulations performed to
                obtain the flux solutions in `fluxes`.

        Returns:
            final_error: mean of the r-squared value of the simulations for the different validation reactions
        """
        for simulation_result, substrate_rate in zip(fluxes, substrate_rates):
            self.parametrization_results.add_fluxes_from_fluxdict(flux_dict=simulation_result,
                                                                  bin_id='final',
                                                                  substrate_uptake_rate= substrate_rate)
        self._init_validation_df(bin_information=[self.min_substrate_uptake_rate, self.max_substrate_uptake_rate])
        validation_results = self.validation_data.sampled_valid_data_df.apply(lambda x: x.abs() if x.dtype!='object' else x)
        error = self._calculate_error_for_reactions(validation_df=validation_results,
                                                    bin_id='final')
        #remove the simulations from the result dataframe
        self.parametrization_results.fluxes_df = self.parametrization_results.fluxes_df[self.parametrization_results.fluxes_df['bin']!='final']
        final_error = np.nanmean(error)
        print('The final error is ', final_error)
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
        esc_in_bin = esc_results_df[esc_results_df['bin'] == bin_id]

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
            if enzyme_id == 'substrate':
                continue

            esc_variability = self._calculate_esc_variability(esc_topn_df, enzyme_id)

            # Determine if ESC variability exceeds the threshold for splitting the bin
            if self._esc_variability_larger_than_threshold(esc_variability):
                self.parametrization_results.bins_to_change.loc[len(self.parametrization_results.bins_to_change)] = [bin_id, True, False]

    def run_genetic_algorithm(self, bin_info: list, esc_topn_df: pd.DataFrame, filename_extension:str) -> None:
        """ Function to run the genetic algorithm.

        Results will be stored in a json, pickle and xlsx file with
        filename defined by: self.hyperparameters.genetic_algorithm_filename_base + filename_extension

        Args:
            bin_info: information about the substrate uptake rates in the bin in the form of:
                        [start, stop, step]
            esc_topn_df: DataFrame containing ESC values, where each column corresponds to an enzyme or substrate
                            uptake reaction, and rows represent different substrate uptake rates
            filename_extension: extension of base filename to save the result of the genetic algorithm
        """
        enzymes_to_evaluate = self._parse_enzymes_to_evaluate(esc_topn_df)
        ga = self._init_genetic_algorithm(bin_info[:2], enzymes_to_evaluate, filename_extension)
        ga.start()

    def restart_genetic_algorithm(self) -> None:
        json_files = self._get_genetic_algorithm_json_files()
        enzymes_to_evaluate = self._determine_enzymes_to_evaluate_for_all_bins(
                                nmbr_kcats_to_pick= self.hyperparameters.number_of_kcats_to_mutate)

        step = (self.max_substrate_uptake_rate - self.min_substrate_uptake_rate) / self.hyperparameters.number_of_bins
        substrate_rates = list(np.arange(self.min_substrate_uptake_rate, self.max_substrate_uptake_rate, step))

        #TODO need to adjust the range of substrate uptake rates for validation
        ga = self._init_genetic_algorithm(substrate_uptake_rates = substrate_rates,
                                          enzymes_to_evaluate = enzymes_to_evaluate,
                                          filename_extension= f'final_run_{self.iteration}')
        ga.restart(json_files)
        return json_files

    def reparametrize_pam(self):
        best_individual_kcat_df, best_indiv_error = self._get_mutated_kcat_values_from_genetic_algorithm()
        error_threshold = self.final_error*(1-self.error_fraction_to_deviate_between_runs)
        if best_indiv_error < error_threshold and best_indiv_error>0:
            return
        for i, row in best_individual_kcat_df.iterrows():
            # current value is coefficient calculated as follows:
            # coeff = 1 / (kcat * 3600 * 1e-6) #3600 to convert to /h to /s *1e-6 to make calculations more accurate
            # need to convert the coefficient to the actual kcat
            # rxn = self.pamodel.reactions.get_by_id(row['rxn_id'])
            # lin_coeff = self.pamodel.constraints[f'EC_{row["id"]}_f'].get_linear_coefficients([rxn.forward_variable])
            # values = [(val, 1/(val*3600*1e-6)) for val in lin_coeff.values()]
            # print(values)
            # self.pamodel.constraints[f'EC_{row["id"]}_f'].set_linear_coefficients({
            #     rxn.forward_variable: (row['value'])})
            new_kcat = 1/(row['value']*3600*1e-6)
            self._change_kcat_value_for_enzyme(enzyme_id= row['id'], kcat_dict = {row['rxn_id']:
                                                                           {'f': new_kcat}
                                                                       })
            # print(self.pamodel.constraints[f'EC_{row["id"]}_f'].get_linear_coefficients([rxn.forward_variable]),
            #       new_kcat, row['value'])

        print('after reparametrization the model is feasible', self._pamodel_is_feasible())


    def save_diagnostics(self, computational_time: float,
                         results_filename:str = None):
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

    def save_final_diagnostics(self, figure:plt.Figure):
        figure.savefig(self.result_figure_file, dpi=100, bbox_inches='tight')
        with pd.ExcelWriter(self.result_diagnostics_file) as writer:
            # Write each DataFrame to a specific sheet
            self.parametrization_results.best_individuals.to_excel(writer, sheet_name='Best_Individuals', index=False)
            self.parametrization_results.computational_time.to_excel(writer, sheet_name='Computational_Time',
                                                                     index=False)
            self.parametrization_results.final_errors.to_excel(writer, sheet_name='Final_Errors', index=False)

    ###########################################################################################################
    #WORKER FUNCTIONS
    ###########################################################################################################
    def _pamodel_is_feasible(self):
        #check for feasibility at a mid substrate uptake rate
        self.pamodel.test((self.max_substrate_uptake_rate -self.min_substrate_uptake_rate/2))
        return self.pamodel.solver.status == 'optimal'

    def _bin_substrate_uptake_rates(self):
        """ Bin the substrate uptake rate in intervals.

         If required, the binsize is adjusted

        Returns:
            dict: dictionary with bin_id: [start, end, stepsize] key value pairs
        """
        substrate_start = self.min_substrate_uptake_rate
        binned_substrate = {}
        bin_range = (self.max_substrate_uptake_rate - self.min_substrate_uptake_rate)/len(self.bins)
        for i in self.bins:
            new_bin = self._make_new_bin(bin_id = i, bin_range = bin_range, substrate_start = substrate_start)
            binned_substrate = {**binned_substrate, **new_bin}
            #update starting concentration for new bin
            substrate_start += bin_range
        return binned_substrate

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
                (bin_id in self.parametrization_results.bins_to_change['bin'].values)):
            new_bin = self._adjust_binsize(bin_id=bin_id, bin_range=bin_range, substrate_start=substrate_start)
        else:
            stepsize = bin_range / self.hyperparameters.bin_resolution
            new_bin = {bin_id: [substrate_start, substrate_start + bin_range, stepsize]}
        return new_bin

    def _adjust_binsize(self, bin_id, bin_range, substrate_start):
        bins_to_change = self.parametrization_results.bins_to_change
        to_split = [(bins_to_change[bins_to_change['bin'] == bin_id]) & (bins_to_change['split'])]
        if len(to_split) > 0:
            stepsize = bin_range * 0.5 / self.hyperparameters.bin_resolution
            return {**{bin_id: [substrate_start, substrate_start + bin_range * 0.5, stepsize]},
                                **{bin_id + 0.1: [substrate_start + bin_range * 0.5, substrate_start + bin_range, stepsize]}}

    def _init_genetic_algorithm(self, substrate_uptake_rates: list,
                                enzymes_to_evaluate: dict,
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
                        'kcat': kcat_value,
                        'sensitivity': esc_value
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
            valid_df = self.validation_data.sampled_valid_data_df,
            filename_save = results_filename,
            substrate_uptake_id = self.substrate_uptake_id, #TODO only get those substrate uptake rate which are included in sampled data
            substrate_uptake_rates = substrate_uptake_rates, # bin_info: [start, stop, step]
            objective_id = self.validation_data._get_biomass_reactions(),
            **self.hyperparameters.genetic_algorithm_hyperparams
        )
        return ga

    def _init_validation_df(self, bin_information):
        validation_results = self.validation_data.valid_data_df
        validation_df = validation_results[
             (abs(validation_results[self.substrate_uptake_id + '_ub']) >= abs(bin_information[0])) &
             (abs(validation_results[self.substrate_uptake_id + '_ub']) <= abs(bin_information[1]))
         ]
        #no data to validate within this range available
        if len(validation_df) == 0: return None
        self.validation_data.sampled_valid_data_df = self._adaptive_sampling(exp_data = validation_df,
                                                        substrate_rxn=self.substrate_uptake_id+'_ub',
                                                        num_samples=self.hyperparameters.bin_resolution,
                                                        min_density=self.hyperparameters.bin_resolution/2,
                                                        max_density=self.hyperparameters.bin_resolution*2)
        # get reference datapoints
        # lower than upper bound and higher than lower bound


    def _init_results_objects(self):
        reactions_to_validate = self.validation_data._reactions_to_validate

        self.parametrization_results.error_df = pd.DataFrame(columns=['bin'] + reactions_to_validate)
        self.parametrization_results.esc_df = pd.DataFrame(columns=['bin', 'substrate', 'enzyme_id', 'rxn_id'])
        self.parametrization_results.fluxes_df = pd.DataFrame(columns=['bin', 'substrate'] + reactions_to_validate)
        self.parametrization_results.initiate_bins_to_change()
        # self.parametrization_results.initiate_result_dfs(reactions_to_validate=reactions_to_validate,
        #                                                   biomass_reaction=self.validation_data._get_biomass_reactions())

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

        esc_results_df['rxn_id'] = esc_results_df['rxn_id'].str.split(',')
        esc_results_df = esc_results_df.explode('rxn_id', ignore_index=True)
        esc_results_df['substrate'] = esc_results_df['substrate'].abs()


        # 1. determine top n values based on average ESC
        # Group by enzyme and calculate the average coefficient for each enzyme
        esc_results_grouped = esc_results_df.groupby('enzyme_id')['coefficient'].agg(['mean']).reset_index()
        esc_results_grouped['absolute_esc_mean'] = esc_results_grouped['mean'].abs()
        # Sort the data based on the average coefficient in descending order
        esc_results_sorted = esc_results_grouped.sort_values(by='absolute_esc_mean', ascending=False)
        # Select the top n enzymes
        top_n_enzymes = esc_results_sorted.head(nmbr_kcats_to_pick)

        # 2. merge with the original DataFrame to connect to the information
        esc_topn_df = pd.merge(top_n_enzymes, esc_results_df.drop_duplicates(), on='enzyme_id')
        # get ids and select those from the dataframe
        esc_topn_df = esc_topn_df.sort_values(by='substrate')
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

        start_esc = esc_df[esc_df['enzyme_id'] == enzyme_id].coefficient.iloc[0]
        end_esc = esc_df[esc_df['enzyme_id'] == enzyme_id].coefficient.iloc[-1]

        # check if the ESCs are nonzero
        if start_esc != 0:
            esc_variability = abs((end_esc - start_esc) / start_esc)
        else:
            # Handle the case when start_esc is zero to avoid division by zero
            esc_variability = 0 if end_esc == 0 else 1  # Assuming that if start_esc is zero, any change is considered 100% change.

        return esc_variability

    def _esc_variability_larger_than_threshold(self, esc_variability: float)-> bool:
        return (esc_variability >= self.hyperparameters.bin_split_deviation_threshold)

    def _calculate_error_for_reactions(self, validation_df: pd.DataFrame,
                                          bin_id: Union[float, int] = None) -> float:
        # calculate error for different exchange rates
        error = []
        for rxn in self.validation_data._reactions_to_validate + self.validation_data._get_biomass_reactions():
            # only select the rows which are filled with data
            validation_data = validation_df.dropna(axis=0, subset=rxn)
            # if there are no reference data points, continue to the next reaction
            if len(validation_data) == 0:
                continue

            flux_df = self.parametrization_results.fluxes_df
            #check if we want to calculate the error for a single bin
            if bin_id is not None: flux_df = flux_df[flux_df['bin'] == bin_id]

            r_squared = self._calculate_r_squared_for_reaction(rxn, validation_df,
                                                               flux_df)
            error += [r_squared]
        return error

    def _calculate_r_squared_for_reaction(self, reaction_id: str, validation_data: pd.DataFrame,
                                          fluxes: pd.DataFrame) -> float:
        substr_rxn = self.substrate_uptake_id+'_ub'
        validation_data[substr_rxn] = [abs(flux) for flux in validation_data[substr_rxn]]
        simulated_data = pd.DataFrame({substr_rxn: [abs(flux) for flux in fluxes['substrate']],
                                       'simulation': fluxes[reaction_id] }) #TODO
        ref_data_rxn = validation_data.merge(simulated_data, on=substr_rxn, how='inner')
        print(simulated_data)
        print(reaction_id)
        print(ref_data_rxn.to_markdown())
        # simulation mean
        data_average = validation_data[reaction_id].mean()
        # error: squared difference
        ref_data_rxn = ref_data_rxn.assign(error=lambda x: (x[reaction_id] - x['simulation']) ** 2)
        print(ref_data_rxn.to_markdown())
        # calculate R^2:
        residual_ss = np.nansum(ref_data_rxn.error)
        total_ss = np.nansum([(data - data_average) ** 2 for data in ref_data_rxn[reaction_id]])
        # calculating r_squared is only feasible of the numerator and the denomenator are both nonzero
        if (residual_ss == 0) | (total_ss == 0):
            r_squared = 0
        else:
            r_squared = 1 - residual_ss / total_ss
        return r_squared

    def _adaptive_sampling(self, exp_data: pd.DataFrame, substrate_rxn: str = 'EX_glc__D_e_ub', num_samples=10, min_density=5,
                          max_density=20):
        """
        Perform adaptive sampling based on data point distribution and variability.

        Parameters:
        - x_exp (array-like): X-coordinates of experimental data.
        - y_exp (array-like): Y-coordinates of experimental data.
        - num_samples (int): Number of samples to generate.
        - min_density (int): Minimum density of samples per unit distance.
        - max_density (int): Maximum density of samples per unit distance.

        Returns:
        - sampled_x (array): Sampled x-coordinates.
        """
        # Calculate distances between consecutive data points
        distances = np.diff(exp_data[substrate_rxn])

        # Compute local data point density as reciprocal of distance
        densities = 1 / distances

        # Normalize densities to [min_density, max_density]
        densities = np.clip(densities, min_density, max_density)

        # Calculate sampling probabilities based on normalized densities
        probabilities = densities / np.sum(densities)

        # Sample indices with replacement based on probabilities
        sampled_indices = np.random.choice(np.arange(len(exp_data) - 1), size=num_samples, p=probabilities)
        exp_data =exp_data.reset_index()
        # Sample x-coordinates from selected indices
        sampled = exp_data.loc[sampled_indices]  # + np.random.uniform(low=0, high=distances[sampled_indices])

        return sampled

    # def _calculate_r_squared_for_reaction(self, reaction_id: str, validation_data: pd.DataFrame,
    #                                       substrate_range: Union[list, np.array],
    #                                       fluxes: pd.DataFrame) -> float:
    #     # calculate difference between simulations and validation data
    #     # get the linear relation of simulation (using the first and the last datapoint)
    #     exp_x = validation_data[self.substrate_uptake_id + '_ub']
    #     exp_y = validation_data[reaction_id]
    #     # Get simulated data points
    #     sim_x = substrate_range
    #     sim_y = fluxes[reaction_id].iloc[0:len(sim_x)]
    #     # Calculate R^2 values using sliding window approach
    #     r_squared = self._interpolated_r_squared(exp_x, exp_y, sim_x, sim_y)
    #     return r_squared


    def _interpolated_r_squared(self, exp_x: list, exp_y: list, sim_x: list, sim_y: list) -> np.array:
        """
        Calculate the R^2 value between experimental and simulated data points within a sliding window.

        Args:
        - exp_x (array-like): X-coordinates of experimental data points.
        - exp_y (array-like): Y-coordinates of experimental data points.
        - sim_x (array-like): X-coordinates of simulated data points.
        - sim_y (array-like): Y-coordinates of simulated data points.
        - window_size (int): Size of the sliding window.

        Returns:
        - r_squared_values (array): squared error values computed for each sliding window.
        """
        error = []
        for i in range(len(exp_x)):
            # Find the indices of experimental data points within the sliding window
            simulation = self._interpolate_datapoint(sim_x, sim_y, exp_x[i])
            error += [(simulation - abs(exp_y[i])) ** 2]
        # calculate R^2:
        exp_average = np.nanmean(exp_y)
        residual_ss = np.nansum(error)
        total_ss = np.nansum([(data - exp_average) ** 2 for data in exp_y])
        # calculating r_squared is only feasible of the numerator and the denomenator are both nonzero
        if (residual_ss == 0) | (total_ss == 0):
            r_squared = 0
        else:
            r_squared = 1 - residual_ss / total_ss

        return r_squared

    def _interpolate_datapoint(self, x_values, y_values, x_target):
        """
        Interpolate a data point between two different x coordinates.

        Args:
        - x_values (array-like): Array of x-coordinates of the data points.
        - y_values (array-like): Array of y-coordinates of the data points.
        - x_target (float): Target x-coordinate to interpolate the data point.

        Returns:
        - interpolated_datapoint (float): Interpolated y-coordinate corresponding to the target x-coordinate.
        """
        # Find the indices of the two nearest x-coordinates
        idx_left = self._find_closest_datapoint(x_values, abs(x_target))
        if x_target < 0 and idx_left > 0:
            idx_right = idx_left - 1
        else:
            idx_right = idx_left + 1
        if idx_right >= len(x_values): idx_right = idx_left - 1

        # Correct datatypes
        if isinstance(x_values, pd.Series): x_values = x_values.values
        if isinstance(y_values, pd.Series): y_values = y_values.values

        # Perform linear interpolation
        x_left, x_right = x_values[idx_left], x_values[idx_right]
        y_left, y_right = y_values[idx_left], y_values[idx_right]
        interpolated_datapoint = y_left + (y_right - y_left) * (abs(x_target) - x_left) / (x_right - x_left)

        return interpolated_datapoint

    def _find_closest_datapoint(self, x_values, target_x):
        """
        Find the data point closest to a specific x value.

        Parameters:
        - x_values (array-like): Array of x-coordinates.
        - target_x (float): Target x value.

        Returns:
        - closest_index (int): Index of the data point closest to the target x value.
        """
        # Calculate absolute differences between target x value and each x-coordinate
        absolute_differences = np.abs(x_values - target_x)
        # Find index of the minimum absolute difference
        closest_index = np.argmin(absolute_differences)
        return closest_index

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
            enzyme_id = row['enzyme_id']
            rxn_id = row['rxn_id']
            rxn = self.pamodel.reactions.get_by_id(rxn_id)
            # print(self.pamodel.constraints[f'EC_{enzyme_id}_f'].get_linear_coefficients([rxn.forward_variable]))
            enzyme_dict = {}
            enzyme_dict['sensitivity'] = row['mean']
            enzyme_dict['reaction'] = rxn_id
            kcat = self.pamodel.constraints[f'EC_{enzyme_id}_f'].get_linear_coefficients([rxn.forward_variable])[rxn.forward_variable]
            # print(self.pamodel.enzymes.get_by_id(enzyme_id).get_kcat_values([rxn_id])['f'])
            enzyme_dict['kcat'] = kcat # TO DO seperate fwd and rev kcat adjustment
            enzymes_to_evaluate[enzyme_id] = enzyme_dict
        return enzymes_to_evaluate

    def _get_genetic_algorithm_json_files(self, subset:str = ''):
        filename_base = self.hyperparameters.genetic_algorithm_filename_base
        folder_path = self.hyperparameters.genetic_algorithm_hyperparams['folderpath_save']

        # Define the regular expression pattern for matching file names
        file_pattern = re.compile(fr'{filename_base}{subset}.*\.json$')
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
        esc_topn_df = self._select_topn_enzymes(esc_results_df,
                                                nmbr_kcats_to_pick)
        enzymes_to_evaluate = self._parse_enzymes_to_evaluate(esc_topn_df)

        return enzymes_to_evaluate

    def _get_mutated_kcat_values_from_genetic_algorithm(self, results_filename:str = None):
        if results_filename is None:
            filename_base = self.hyperparameters.genetic_algorithm_filename_base
            folder_path = self.hyperparameters.genetic_algorithm_hyperparams['folderpath_save']

            results_filename = os.path.join(folder_path,
                                                 filename_base + f'final_run_{self.iteration}.xlsx')
        final_run_results_best_indiv = pd.read_excel(results_filename, sheet_name='best_individual').drop(['type'],
                                                                                                          axis = 1)
        error_df = pd.read_excel(results_filename, sheet_name= 'final_population')
        error = float(error_df.iloc[0]['fitness_weighted_sum'])
        return final_run_results_best_indiv ,error

    def _change_kcat_value_for_enzyme(self, enzyme_id:str, kcat_dict:dict) -> None:
        self.pamodel.change_kcat_value(enzyme_id=enzyme_id, kcats=kcat_dict)

    def _remove_result_files(self, file_base: Union[list, str]) -> None:
        """ Removes files resulting from genetic algorithm runs.

        Args:
            file_base (list, string): list of stings or pathlib.Path objects or a string with the files
            to remove. If there are file extensions, they will be removed.
        """
        if not hasattr(file_base, '__iter__'): file_base = [file_base]
        for file in file_base:
            if not isinstance(file, str): file = str(file)
            file_path_base = file.split('.')[0]
            [os.remove(file_path_base + file_type) for file_type in ['.json', '.xlsx', '.pickle']]

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
        enz_rxn_kcat_row['value'] = 1/(enz_rxn_kcat_row['value']*3600*1e-6)
        return enz_rxn_kcat_row.to_list()


    #################################################################################################################
    # VISUALIZATION
    #################################################################################################################

    def plot_valid_data(self):
        # plot flux changes with glucose uptake
        fig, axs = plt.subplots(2,2, dpi=100)

        for r, ax in zip(self.validation_data._reactions_to_plot, axs.flatten()):
            # plot data
            x = [abs(glc) for glc in self.validation_data.valid_data_df[self.substrate_uptake_id +'_ub']]
            y = [abs(data) for data in self.validation_data.valid_data_df[r]]
            ax.set_ylabel(r)
            ax.scatter(x, y,
                           color='black', marker='o', s=30, linewidths=1.3,
                           facecolors=None, zorder=0,
                           label='Data')
            ax.set_xlabel(self.substrate_uptake_id + ' $mmol/g_{CDW}/h$')
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
            cmap = plt.get_cmap('viridis')
            color = to_hex(cmap(self.iteration/(self.hyperparameters.threshold_iteration-1)))

        fluxes, substrate_range = self.run_simulations_to_plot(save_esc)

        for r, ax in zip(self.validation_data._reactions_to_plot, axs.flatten()):
            # plot data
            line = ax.plot(substrate_range, [abs(f[r]) for f in fluxes], linewidth=2.5,
                           zorder=5, color=color)

        fig.canvas.draw()
        fig.canvas.flush_events()
        if return_fluxes: return fig, fluxes, substrate_range
        return fig

    def run_simulations_to_plot(self, substrate_rates: Union[np.array, list, pd.Series] = None,
                                save_fluxes_esc:bool = False):
        fluxes = list()
        substrate_range = list()
        if substrate_rates is None:
            step = (self.max_substrate_uptake_rate-self.min_substrate_uptake_rate)/10
            substrate_rates = np.arange(self.min_substrate_uptake_rate, self.max_substrate_uptake_rate, step)
        for substrate in substrate_rates:
            if substrate>=0:
                self.pamodel.change_reaction_bounds(rxn_id=self.substrate_uptake_id,
                                                lower_bound=0, upper_bound=substrate)
            else:
                self.pamodel.change_reaction_bounds(rxn_id=self.substrate_uptake_id,
                                                    lower_bound=substrate, upper_bound=0)
            # solve the model
            sol_pam = self.pamodel.optimize()

            if self.pamodel.solver.status == 'optimal' and self.pamodel.objective.value != 0:
                substrate_range += [abs(substrate)]
                fluxes.append(sol_pam.fluxes)
                if save_fluxes_esc: self.save_pamodel_simulation_results(substrate_uptake_rate=substrate,
                                                                         bin_id= 'no bins')
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
        clustering = AgglomerativeClustering(n_clusters=n_clusters, linkage='ward')

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