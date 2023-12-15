from typing import Union
import numpy as np
import os
import time
import matplotlib.pyplot as plt
from PAModelpy.PAModel import PAModel

from .PAM_data_classes import ValidationData, HyperParameters, ParametrizationResults

class PAMParametrizer():
    def __init__(self, pamodel:PAModel,
                 validation_data: ValidationData = ValidationData,
                 hyperparameters: HyperParameters = HyperParameters,
                 max_substrate_uptake_rate: Union[float, int] = 11,
                 min_substrate_uptake_rate: Union[float, int] = 0):
        self.core_genetic_algorithm = None
        self.pamodel = pamodel
        self.validation_data = validation_data
        self.hyperparameters = hyperparameters
        self.parametrization_results = ParametrizationResults
        self.enzyme_ids = [enzyme.id for enzyme in pamodel.enzymes]

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
            # 0. initiate structures to save results and keep track of process
            self.iteration +=1
            self.parametrization_results.initiate_result_dfs(enzyme_ids=self.enzyme_ids,
                                                             reactions_to_validate= self.validation_data.reactions_to_validate,
                                                             growth_rate= self.validation_data.growth_rates)

            # 1. adjust binsize if required and get the binned substrate values
            binned_substrate = self.bin_substrate_uptake_rates()


            # 2. Run model in bins, get sensitivities and calculate errors
            for index, bin in binned_substrate.items():
                print(
                    f'The following range of substrate uptake rates will be analyzed: {bin[0]} - {bin[1]} mmol/g_cdw/h, with steps of {bin[2]} mmol/g_cdw/h')
                for substrate in np.arange(bin[0], bin[1], bin[2]):
                    self.run_simulations(substrate, index)
                # calculate the error for the different exchange rates
                self.calculate_error(bin, index, substrate)
                #print running time to check on progress
                print('time elapsed: ', time.perf_counter() - start, 'sec, ', (time.perf_counter() - start) / 60, 'min',
                  (time.perf_counter() - start) / 3600, 'hour\n')

            #get the most sensitive enzymes and if an enzyme is in the topn enzymes in multiple bins, sample 1 error/sensitivity
            self.compare_sensitivities()
            self.sample_duplicate_enzymes()
            # 6. Adjust sensitive enzymes using: i. E/kcat*sens, ii. Genetic algorithm like with a probability distribution of picking, iii. reinforcement learning
            self.adjust_kcats_gradientdescent()

            #visualize results
            fig = self.plot_simulation(fig, axs)

            # 6. Display progress and repeat
            print('time elapsed: ', time.perf_counter() - start, 'sec, ', (time.perf_counter() - start) / 60, 'min',
                  (time.perf_counter() - start) / 3600, 'hour')
            print('Done with iteration number ', self.iteration)
            print('-------------------------------------------------------------------------------------------')

        fig.savefig(self.result_figure_file, dpi=100, bbox_inches='tight')

    def run_genetic_algorithm(self):
        pass

    def calculate_error(self):
        pass

    def determine_most_sensitive_enzymes(self):
        pass

    def reparametrize(self):
        pass

    def run_pamodel(self):
        pass

    def save_diagnostics(self):
        pass

    ###########################################################################################################
    #WORKER FUNCTIONS
    ###########################################################################################################

    def bin_substrate_uptake_rates(self):
        """
        Bin the substrate uptake rate in intervals. If required, the binsize is adjusted
        :return: dict: dictionary with bin_id: [start, end, stepsize] key value pairs
        """
        substrate_start = self.min_substrate_uptake_rate
        binned_substrate = {}
        bin_range = (self.max_substrate_uptake_rate - self.min_substrate_uptake_rate)/len(self.bins)
        for i in self.bins:
            new_bin = self.make_new_bin(bin_id = i, bin_range = bin_range, substrate_start = substrate_start)
            binned_substrate = {**binned_substrate, **new_bin}
            #update starting concentration for new bin
            substrate_start += bin_range
        return binned_substrate

    def make_new_bin(self, bin_id:Union[float, int], bin_range: Union[float, int],
                     substrate_start:Union[float, int]):
        """
        Creating a dictionary with all information about a single bin (range of substrate uptake rates).
        Per bin it saves the start, end and stepsize.
        The function also checks if the bin should be adapted based on the previous iteration of the workflow

        :param bin_id: identifier of the bin to make
        :param bin_range: length of the bin in the units of the substrate uptake rate
        :param substrate_start: the first substrate uptake rate in the bin
        :return: dict: dictionary with bin_id: [start, end, stepsize] key value pairs
        """
        # if bins should be smaller, adjust binsize
        if ((len(self.parametrization_results.bins_to_change) > 0) and
                (bin_id in self.parametrization_results.bins_to_change['bin'].values)):
            new_bin = self.adjust_binsize(bin_id=bin_id, bin_range=bin_range, substrate_start=substrate_start)
        else:
            stepsize = bin_range / self.hyperparameters.bin_resolution
            new_bin = {bin_id: [substrate_start, substrate_start + bin_range, stepsize]}
        return new_bin

    def adjust_binsize(self, bin_id, bin_range, substrate_start):
        bins_to_change = self.parametrization_results.bins_to_change
        to_split = [(bins_to_change[bins_to_change['bin'] == bin_id]) & (bins_to_change['split'])]
        if len(to_split) > 0:
            stepsize = bin_range * 0.5 / self.hyperparameters.bin_resolution
            return {**{bin_id: [substrate_start, substrate_start + bin_range * 0.5, stepsize]},
                                **{bin_id + 0.1: [substrate_start + bin_range * 0.5, substrate_start + bin_range, stepsize]}}
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

        for r, ax in zip(self.validation_data.reactions_to_plot, axs.flatten()):
            # plot data
            line = ax.plot(self.substrate_range, [abs(f[r]) for f in self.fluxes], linewidth=2.5,
                           zorder=5, color=f'#{self.color}')

        fig.canvas.draw()
        fig.canvas.flush_events()
        return fig

