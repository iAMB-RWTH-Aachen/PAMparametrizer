from typing import Union
import numpy as np
import os
import time
import matplotlib.pyplot as plt
from PAModelpy.PAModel import PAModel
import pandas as pd
from scipy.stats import linregress
import re

from .PAM_data_classes import ValidationData, HyperParameters, ParametrizationResults
from Modules.genetic_algorithm_parametrization import GAPOGaussian as GAPOGauss


class PAMParametrizer():
    def __init__(self, pamodel:PAModel,
                 validation_data: ValidationData = ValidationData,
                 hyperparameters: HyperParameters = HyperParameters,
                 substrate_uptake_id: str = 'EX_glc__D_e',
                 max_substrate_uptake_rate: Union[float, int] = 11,
                 min_substrate_uptake_rate: Union[float, int] = 0):

        self.core_genetic_algorithm = None
        self.pamodel = pamodel
        self.validation_data = validation_data
        self.hyperparameters = hyperparameters
        self.parametrization_results = ParametrizationResults()
        self.enzyme_ids = [enzyme.id for enzyme in pamodel.enzymes]

        self.substrate_uptake_id = substrate_uptake_id
        self.min_substrate_uptake_rate = min_substrate_uptake_rate
        self.max_substrate_uptake_rate = max_substrate_uptake_rate

        # attributes for keeping track of the workflow
        self.iteration =0
        self.bins = list(range(self.hyperparameters.number_of_bins))
        self.result_figure_file = os.path.join(os.path.split(os.path.split(os.getcwd())[0])[0], 'Results', 'pam_parametrizer.png')


    def run(self):
        #setup plot to visualize progress
        fig, axs = self.plot_valid_data()

        #keep track of time for computational performance
        start = time.perf_counter()

        #perform parametrization until a max number of iterations
        while self.iteration <= self.hyperparameters.threshold_iteration:
            # 0. initiate structures and objects to save results and keep track of process
            self.iteration +=1
            self._init_results_objects()

            # 1. Get the binned substrate values and adjust binsize if required
            binned_substrate = self._bin_substrate_uptake_rates()

            # 2. Run model in bins, get sensitivities and calculate errors
            for bin_id, bin in binned_substrate.items():
                print(
                    f'The following range of substrate uptake rates will be analyzed: {bin[0]} - {bin[1]} '
                    f'mmol/g_cdw/h, with steps of {bin[2]} mmol/g_cdw/h')
                self.run_pamodel_simulations_in_bin(bin_information=bin)
                # calculate the error for the different exchange rates
                self.calculate_error(bin, bin_id)
                # if there is no solution, continue to the next bin
                if len(self.parametrization_results.esc_df[self.parametrization_results.esc_df['bin'] == bin_id]) == 0:
                    continue
                # get the most sensitive enzymes
                top_enzyme_sensitivities = self.determine_most_sensitive_enzymes(bin_id,
                                                                                 self.hyperparameters.number_of_kcats_to_mutate)
                self.determine_bin_to_split(top_enzyme_sensitivities, bin_id)

                self.run_genetic_algorithm(bin_info = bin,
                                           esc_topn_df= top_enzyme_sensitivities,
                                           filename_extension = str(bin_id))

                #print running time to check on progress
                print('time elapsed: ', time.perf_counter() - start, 'sec, ', (time.perf_counter() - start) / 60, 'min',
                  (time.perf_counter() - start) / 3600, 'hour\n')

            self.restart_genetic_algorithm()
            self.reparametrize()

            #visualize results
            fig = self.plot_simulation(fig, axs)

            # 6. Display progress and repeat
            print('time elapsed: ', time.perf_counter() - start, 'sec, ', (time.perf_counter() - start) / 60, 'min',
                  (time.perf_counter() - start) / 3600, 'hour')
            print('Done with iteration number ', self.iteration)
            print('-------------------------------------------------------------------------------------------')

        fig.savefig(self.result_figure_file, dpi=100, bbox_inches='tight')

    def run_pamodel_simulations_in_bin(self, bin_information:dict) -> None:
        """
        Use the range of substrate uptake rate indicated in the bin_information to run simulations with the PAModel
        :param bin_information: dictionary with bin_id:[start, stop, step] key:value pairs, where start, stop and step
                            relate to the start, end and stepsize of the substrate uptake rate range
        """
        # self._init_results_objects()
        for bin_id, bin_info in bin_information.items():
            start, stop, step = bin_info[0], bin_info[1], bin_info[2]
            print(
                f'The following range of substrate uptake rates will be analyzed: {start - stop} '
                f'mmol/g_cdw/h, with steps of {step} mmol/g_cdw/h')
            for substrate_uptake_rate in np.arange(start, stop, step):
                self.run_pamodel_simulation(substrate_uptake_rate, bin_id)

    def run_pamodel_simulation(self, substrate_uptake_rate: Union[float, int],
                               bin_id: Union[str, float, int]) -> None:
        """
        Running PAModel simulations and saving the resulting fluxes and enzymes sensitivities coefficient
        for later analysis

        :param substrate_uptake_rate: used to constrain the PAModel
        :param bin_id: identifier of the bin in which this simulation is run (for saving purposes)
        """

        print('Substrate uptake rate ', substrate_uptake_rate, ' mmol/gcdw/h')
        with self.pamodel:
            # change glucose uptake rate
            self.pamodel.change_reaction_bounds(rxn_id=self.substrate_uptake_id,
                                                lower_bound=0, upper_bound=substrate_uptake_rate)
            # solve the model
            self.pamodel.optimize()

            if self.pamodel.solver.status == 'optimal' and self.pamodel.objective.value != 0:
                self.save_pamodel_simulation_results(substrate_uptake_rate, bin_id)

    def save_pamodel_simulation_results(self, substrate_uptake_rate: Union[float, int],
                                        bin_id: Union[str, float, int]) -> None:
        """
            Saving the resulting fluxes and enzymes sensitivities coefficient from a successful PAModel simulation
            for later analysis

            :param substrate_uptake_rate: used to constrain the PAModel
            :param bin_id: identifier of the bin in which this simulation is run (for saving purposes)
        """

        self.parametrization_results.substrate_range += [substrate_uptake_rate]
        self.parametrization_results.add_fluxes(self.pamodel, bin_id, substrate_uptake_rate)
        self.parametrization_results.add_enzyme_sensitivity_coefficients(self.pamodel.enzyme_sensitivity_coefficients,
                                                                          bin_id, substrate_uptake_rate)

    def calculate_error(self, bin_information: list, bin_id: Union[float, int]) -> None:
        """
        Evaluate the model simulations compared to a reference dataset within a specified substrate uptake range.

        This method calculates the average difference between model simulations and reference data for the total substrate
        uptake range and available reactions within the specified bin.

        Notes:
        - The reference dataset is filtered to include only data points within the specified substrate uptake range.
        - Errors are calculated for different exchange rates and biomass reactions.
        - The '_calculate_r_squared_for_reaction' method is used to determine the R-squared value for each reaction.
        - The calculated errors are stored in the 'error_df' attribute for further analysis.

        :param: bin_information: List containing upper and lower bounds of the substrate uptake range for the bin.
        :param: bin_id: Identifier for the bin under consideration.
        """

        validation_results = self.validation_data.valid_data_df
        error = []
        # get reference datapoints
        # lower than upper bound and higher than lower bound
        validation_df = validation_results[
             (validation_results[self.substrate_uptake_id] <= -bin_information[0]) &
             (validation_results[self.substrate_uptake_id] >= -bin_information[1])
         ]

        # calculate error for different exchange rates
        for rxn in self.validation_data._reactions_to_validate + self.validation_data._get_biomass_reactions():
            # only select the rows which are filled with data
            validation_data = validation_df.dropna(axis=0, subset=rxn)
            # if there are no reference data points, continue to the next reaction
            if len(validation_data) == 0:
                continue

            r_squared = self._calculate_r_squared_for_reaction(rxn, validation_df, bin_information, bin_id)

            error += [r_squared]

        error_df = self.parametrization_results.error_df
        self.parametrization_results.error_df.loc[len(error_df)] = [bin_id] + error

    def determine_most_sensitive_enzymes(self, bin_id: Union[float, int], nmbr_kcats_to_pick: int) -> pd.DataFrame:
        """
        Determine the top n (with n being the nmbr_of_kcats_to_pick) sensitive enzymes in a specific
        bin and save them to the dataclasses stored in parametrization_results. Selection is
        based on the enzyme sensitivity coefficients (refer to PAModelpy documentation for more
        details on calculations)

        :param bin_id: bin identifier
        :param nmbr_kcats_to_pick: number of enzymes to select for optimization

        :return esc_topn_df: DataFrame containing the (mean) ESC values of the selected enzymes, with as columns:
                            [enzyme_id,  mean,  absolute_esc_mean,  bin,  substrate, rxn_id,  coefficient]
        """
        esc_results_df = self.parametrization_results.esc_df
        # 0. get enzyme sensitivities
        esc_in_bin = esc_results_df[esc_results_df['bin'] == bin_id]

        # 1. determine top n values based on average ESC
        # Group by enzyme and calculate the average coefficient for each enzyme
        esc_in_bin_grouped = esc_in_bin.groupby('enzyme_id')['coefficient'].agg(['mean']).reset_index()
        esc_in_bin_grouped['absolute_esc_mean'] = esc_in_bin_grouped['mean'].abs()
        # Sort the data based on the average coefficient in descending order
        esc_in_bin_sorted = esc_in_bin_grouped.sort_values(by='absolute_esc_mean', ascending=False)
        # Select the top n enzymes
        top_n_enzymes = esc_in_bin_sorted.head(nmbr_kcats_to_pick)

        # 2. merge with the original DataFrame to connect to the information
        esc_topn_df = pd.merge(top_n_enzymes, esc_in_bin.drop_duplicates(), on='enzyme_id')
        # get ids and select those from the dataframe
        # topn_ids = esc_topn_df.enzyme_id.to_list()
        esc_topn_df = esc_topn_df.sort_values(by='substrate')

        # for id in topn_ids:
        #     mean_esc = esc_topn_df[esc_topn_df['enzyme_id' == id]]['mean']
        #     esc_results_df.loc[len(esc_results_df)] = [bin_id, mean_esc, id]
        return esc_topn_df

    def determine_bin_to_split(self, esc_topn_df: pd.DataFrame, bin_id: Union[float, int]) -> None:
        """
        Determine whether a bin should be split based on the variability of ESC values.

        This function calculates the variability of ESC values within a bin and checks if it exceeds a certain threshold,
        indicating a significant change in the enzymatic substrate concentrations (ESCs) over the course of the bin.

        :param esc_topn_df: DataFrame containing ESC values, where each column corresponds to an enzyme or substrate
            uptake rate, and rows represent different substrate uptake rates
        :param bin_id: Identifier for the bin under consideration.
        """

        for enzyme_id in esc_topn_df.columns:
            if enzyme_id == 'substrate':
                continue

            esc_variability = self._calculate_esc_variability(esc_topn_df, enzyme_id)

            # Determine if ESC variability exceeds the threshold for splitting the bin
            if self._esc_variability_larger_than_threshold(esc_variability):
                self.parametrization_results.bins_to_change.loc[len(self.parametrization_results.bins_to_change)] = [bin_id, True, False]

    def run_genetic_algorithm(self, bin_info:list, esc_topn_df: pd.DataFrame, filename_extension:str) -> None:
        """
        Function to run the genetic algorithm. Results will be stored in a json, pickle and xlsx file with
        filename defined by: self.hyperparameters.genetic_algorithm_filename_base + filename_extension

        :param bin_info: information about the substrate uptake rates in the bin in the form of:
                        [start, stop, step]
        :param esc_topn_df: DataFrame containing ESC values, where each column corresponds to an enzyme or substrate
                            uptake reaction, and rows represent different substrate uptake rates
        :param filename_extension: extension of base filename to save the result of the genetic algorithm
        """
        enzymes_to_evaluate = self._parse_enzymes_to_evaluate(esc_topn_df)

        ga = self._init_genetic_algorithm(bin_info, enzymes_to_evaluate, filename_extension)
        ga.start()

    def restart_genetic_algorithm(self):
        json_files = self._get_genetic_algorithm_json_files()
        enzymes_to_evaluate  = self._determine_enzymes_to_evaluate_from_all_bins()
        substrate_start = self.validation_data.valid_data_df[self.substrate_uptake_id].iloc[0]
        substrate_stop = self.validation_data.valid_data_df[self.substrate_uptake_id].iloc[-1]

        ga = self._init_genetic_algorithm(substrate_uptake_rates = [substrate_start, substrate_stop],
                                          enzymes_to_evaluate = enzymes_to_evaluate,
                                          filename_extension= f'_final_run_{self.iteration}')
        ga.restart(json_files)


    def reparametrize(self):
        pass

    def save_diagnostics(self):
        pass

    ###########################################################################################################
    #WORKER FUNCTIONS
    ###########################################################################################################

    def _bin_substrate_uptake_rates(self):
        """
        Bin the substrate uptake rate in intervals. If required, the binsize is adjusted
        :return: dict: dictionary with bin_id: [start, end, stepsize] key value pairs
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
        """
        Creating a dictionary with all information about a single bin (range of substrate uptake rates).
        Saves the start, end and stepsize per bin.
        The function also checks if the bin should be adapted based on the previous iteration of the workflow

        :param bin_id: identifier of the bin to make
        :param bin_range: length of the bin in the units of the substrate uptake rate
        :param substrate_start: the first substrate uptake rate in the bin
        :return: dict: dictionary with bin_id: [start, end, stepsize] key value pairs
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
        Initiate core genetic algorithm object. Upon initiation the core genetic algorithm will
        initiate the genetic algorithm runner (ga_param) and the fitness evaluation class.
        Hyperparameters related of the genetic algorithm can be adjusted in
        self.hyperparameters.genetic_algorithm_hyperparameters (dict). The sampling method (uniform
        or from a normal(Gaussian) distribution) can be chosen by changing the core genetic algorithm
        object in self.hyperparameters.genetic_algorithm.

        Files will be saved with filename defined by:
        self.hyperparameters.genetic_algorithm_filename_base + filename_extension

        :param substrate_uptake_rates: information about the substrate uptake rates to evaluate in the form of:
                        [start, stop]
        :param enzymes_to_evaluate: dictionary with required information about which enzyme parameters
                                    should be optimized. Format:
                                    {'enzyme_id':{
                                        'reaction':'reaction_id',
                                        'kcat':kcat_value,
                                        'sensitivity': esc_value
                                                }}
        :param filename_extension: extension of base filename to save the result of the genetic algorithm
        :return: ga: core genetic algorithm object ready to be started
        """

        results_filename = self.hyperparameters.genetic_algorithm_filename_base + filename_extension

        ga = self.hyperparameters.genetic_algorithm(
            model=self.pamodel,
            enzymes_to_eval= enzymes_to_evaluate,
            valid_df = self.validation_data.valid_data_df,
            filename_save = results_filename,
            substrate_uptake_id = self.substrate_uptake_id,
            substrate_uptake_rates = [substrate_uptake_rates[0], substrate_uptake_rates[1]], # bin_info: [start, stop, step]
            objective_id = self.validation_data._get_biomass_reactions(),
            **self.hyperparameters.genetic_algorithm_hyperparams
        )
        return ga

    def _init_results_objects(self):
        self.parametrization_results.initiate_result_dfs(reactions_to_validate=self.validation_data._get_reactions_to_validate(),
                                                         biomass_reaction=self.validation_data._get_biomass_reactions())

    def _calculate_esc_variability(self, esc_df:pd.DataFrame ,enzyme_id: str) -> float:
        """
        ESC variability is described as a fractional change.
        :param esc_df: DataFrame containing ESC values.
        :param enzyme_id: ID of the enzyme.
        :return: ESC variability as a percentage change.
        """
        start_esc = esc_df[enzyme_id].iloc[0]
        end_esc = esc_df[enzyme_id].iloc[-1]

        # check if the ESCs are nonzero
        if start_esc != 0:
            esc_variability = abs((end_esc - start_esc) / start_esc)
        else:
            # Handle the case when start_esc is zero to avoid division by zero
            esc_variability = 0 if end_esc == 0 else 1  # Assuming that if start_esc is zero, any change is considered 100% change.

        return esc_variability

    def _esc_variability_larger_than_threshold(self, esc_variability: float)-> bool:
        return (esc_variability >= self.hyperparameters.bin_split_deviation_threshold)

    def _calculate_r_squared_for_reaction(self, reaction_id: str, validation_data: pd.DataFrame, bin_information: list,
                                          bin_id: Union[float, int]) -> float:
        flux_df = self.parametrization_results.fluxes_df
        fluxes = flux_df[flux_df['bin'] == bin_id]

        # calculate difference between simulations and validation data
        # get the linear relation of simulation (using the first and the last datapoint)
        line = linregress(x=[abs(bin_information[0]), abs(bin_information[1])], y=[fluxes[reaction_id].iloc[0], fluxes[reaction_id].iloc[-1]])
        ref_data_rxn = validation_data.assign(
                simulation=lambda x: line.intercept + line.slope * x[self.substrate_uptake_id])

        # simulation mean
        data_average = ref_data_rxn[reaction_id].mean()
        # error: squared difference
        ref_data_rxn = ref_data_rxn.assign(error=lambda x: (x[reaction_id] - x['simulation']) ** 2)

        # calculate R^2:
        residual_ss = sum(ref_data_rxn.error)
        total_ss = sum([(data - data_average) ** 2 for data in ref_data_rxn[reaction_id]])
        r_squared = 1 - residual_ss / total_ss
        return r_squared

    def _parse_enzymes_to_evaluate(self, esc_topn_df: pd.DataFrame) -> dict:
        """
        Function to parse the enzyme-reaction pairs for which the kinetic parameters should be adjusted
        from the enzyme sensitivity coefficients outputted by the model to an input of the genetic algorithm.

        :param esc_topn_df: DataFrame containing the (mean) ESC values of the selected enzymes, with as columns:
                            [enzyme_id,  mean,  absolute_esc_mean,  bin,  substrate, rxn_id,  coefficient]
        :return: enzymes_to_evaluate: dictionary with required information about which enzyme parameters
                                    should be optimized. Format:
                                    {'enzyme_id':{
                                        'reaction':'reaction_id',
                                        'kcat':kcat_value,
                                        'sensitivity': esc_value
                                                }}
        """
        enzymes_to_evaluate = {}

        for index, row in esc_topn_df.iterrows():
            enzyme_id = row['enzyme_id']
            rxn_id = row['rxn_id']
            enzyme_dict = {}
            enzyme_dict['sensitivity'] = row['mean']
            enzyme_dict['reaction'] = rxn_id
            enzyme_dict['kcat'] = self.pamodel.enzymes.get_by_id(enzyme_id).get_kcat_values([rxn_id])['f'] # TO DO seperate fwd and rev kcat adjustment
            enzymes_to_evaluate[enzyme_id] = enzyme_dict
        return enzymes_to_evaluate

    def _get_genetic_algorithm_json_files(self):
        filename_base = self.hyperparameters.genetic_algorithm_filename_base
        folder_path = self.hyperparameters.genetic_algorithm_hyperparams['folderpath_save']

        # Define the regular expression pattern for matching file names
        file_pattern = re.compile(fr'{filename_base}.*\.json$')

        # Get a list of all files in the folder and filter using the regular expression
        json_files = [f for f in os.listdir(folder_path) if file_pattern.match(f)]
        return json_files

    def _determine_enzymes_to_evaluate_from_all_bins(self) -> dict:
        #TODO
        return dict

    def plot_valid_data(self):
        # plot flux changes with glucose uptake
        fig, axs = plt.subplots(2,2, dpi=100)

        for r, ax in zip(self.validation_data.reactions_to_validate, axs.flatten()):
            # plot data
            x = [-glc for glc in self.validation_data.valid_data_df[self.substrate_uptake_id]]
            y = [abs(data) for data in self.validation_data.valid_data_df[r]]
            ax.set_ylabel(r)
            ax.scatter(x, y,
                           color='firebrick', marker='o', s=30, linewidths=1.3,
                           facecolors=None, zorder=0,
                           label='Data')
        plt.ion()   # set interactive mode
        fig.tight_layout()
        fig.show()

        return fig, axs

    def plot_simulation(self, fig, axs):
        #adjust color to visualize progress
        self.color -= 50

        for r, ax in zip(self.validation_data._reactions_to_plot, axs.flatten()):
            # plot data
            line = ax.plot(self.parametrization_results.substrate_range, [abs(f[r]) for f in self.fluxes], linewidth=2.5,
                           zorder=5, color=f'#{self.color}')

        fig.canvas.draw()
        fig.canvas.flush_events()
        return fig
