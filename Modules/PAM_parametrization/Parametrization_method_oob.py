#package
import os
import numpy as np
import pandas as pd
from scipy.stats import linregress
import sys
import time
import matplotlib.pyplot as plt
import tabulate

#parallel processing
import multiprocessing as mp

from PAModelpy.configuration import Config
sys.path.append('../Examples/')
from example_PAM_generation_newparams import set_up_pam
#also need xlwt for to_excel function in pandas

rxn_id = [Config.ACETATE_EXCRETION_RXNID, Config.CO2_EXHANGE_RXNID, Config.OXYGEN_UPTAKE_RXNID, Config.BIOMASS_REACTION]


class Parametrizer:

    #constants
    CPU = mp.cpu_count()-2
    SUBSTRATE_END = 11 #mmol/gcdw/h
    NMBR_BINS = 3
    BIN_RESOLUTION = 5
    NMBR_KCAT = 10 #to change
    CHANGE_THRES = 0.20 #faction of allowed change within a bin
    GLUCOSE_EXCHANGE_RXNID = Config.GLUCOSE_EXCHANGE_RXNID
    LEARNING_RATE = 10
    START_COLOR = 440154

    def __init__(self, valid_data_file: pd.DataFrame, valid_rxns: list,
                 substrate_uptake_id:str = 'EX_glc__D_e',iterations:int = 10):
        self.valid_data_df = valid_data_file
        self.valid_rxn_ids = valid_rxns

        #only get exchanges and growth rate
        self.growth_rate = [data for data in valid_data_file.columns if data.split('_')[0]=='BIOMASS']
        self.reactions_with_data = [data for data in valid_data_file.columns if data.split('_')[0]=='EX']

        #the model
        self.pamodel =set_up_pam()

        #information for parametrization loop
        self.max_iterations = iterations
        self.substrate_uptake_id = substrate_uptake_id
        self.max_substrate_uptake = self.SUBSTRATE_END
        self.nmbr_bins = self.NMBR_BINS
        self.bin_resolution = self.BIN_RESOLUTION
        self.bins = list(range(self.nmbr_bins))
        self.bins_to_change = pd.DataFrame(columns=['bin', 'split', 'merge'])
        self.change_threshold = self.CHANGE_THRES
        self.nmbr_kcats_to_change = self.NMBR_KCAT

        #objects for storing results
        self.error_df = pd.DataFrame(columns=['bin', 'substrate'] + self.reactions_with_data + self.growth_rate)
        self.enzyme_ids = [enz.id for enz in list(self.pamodel.enzymes)]
        self.fac_df = pd.DataFrame(columns=['bin', 'substrate'] + self.enzyme_ids)
        self.sensitive_enzymes = pd.DataFrame(columns=['bin', 'mean_sensitivity', 'enzyme_id'])
        self.fluxes = []
        self.substrate_range = []

        #for plotting
        self.color = self.START_COLOR

        self.iteration = 0

    def plot_valid_data(self):
        # plot flux changes with glucose uptake
        fig, axs = plt.subplots(2,2, dpi=100)

        for r, ax in zip(self.valid_rxn_ids, axs.flatten()):
            # plot data
            x = [-glc for glc in self.valid_data_df[self.substrate_uptake_id]]
            y = [abs(data) for data in self.valid_data_df[r]]
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

        for r, ax in zip(rxn_id, axs.flatten()):
            # plot data
            line = ax.plot(self.substrate_range, [abs(f[r]) for f in self.fluxes], linewidth=2.5,
                           zorder=5, color=f'#{self.color}')

        fig.canvas.draw()
        fig.canvas.flush_events()
        return fig

    def run(self):
        #setup plot to visualize progress
        fig, axs = self.plot_valid_data()

        start=time.perf_counter()
        #1. set up model with dummy parameters
        pamodel = set_up_pam()
        enzyme_ids = [enz.id for enz in list(pamodel.enzymes)]
        while self.iteration <= self.max_iterations:
            #keep count of iteration
            self.iteration +=1
            #2. adjust binsize if required and get the binned substrate values
            binned_substrate = self.adjust_binsize()

            #empty result storage dataframes
            self.error_df = pd.DataFrame(columns=['bin', 'substrate'] + self.reactions_with_data + self.growth_rate)
            self.fac_df = pd.DataFrame(columns=['bin', 'substrate'] + enzyme_ids)

            # 3. Run model in bins, get sensitivities and calculate errors
            for index, bin in binned_substrate.items():
                print(
                    f'The following range of substrate uptake rates will be analyzed: {bin[0]} - {bin[1]} mmol/g_cdw/h, with steps of {bin[2]} mmol/g_cdw/h')
                for substrate in np.arange(bin[0], bin[1], bin[2]):
                    self.run_simulations(substrate, index)
                # calculate the error for the different exchange rates
                self.calculate_r_squared(bin, index, substrate)
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

        fig.savefig('parametrizationfig.png', dpi=100, bbox_inches='tight')

    def adjust_binsize(self):
        substrate_start = 0
        binned_substrate = {}
        bin_diff = self.max_substrate_uptake/len(self.bins)
        for i in self.bins:
            # if bins should be smaller, adjust binsize
            if len(self.bins_to_change) > 0 and i in self.bins_to_change['bin']:
                to_split = self.bins_to_change[(self.bins_to_change['bin'] == i) & (self.bins_to_change['split'])]
                if len(to_split) > 0:
                    stepsize = bin_diff * 0.5 / self.bin_resolution
                    binned_substrate = {**binned_substrate, **{i: [substrate_start, substrate_start + bin_diff * 0.5, stepsize]},
                                  **{i + 0.1: [substrate_start + bin_diff * 0.5, substrate_start + bin_diff, stepsize]}}
            else:
                stepsize = bin_diff / self.bin_resolution
                binned_substrate = {**binned_substrate, **{i: [substrate_start, substrate_start + bin_diff, stepsize]}}
            substrate_start += bin_diff
        return binned_substrate

    def run_simulations(self, substrate, index):
        print('Substrate uptake rate ', substrate, ' mmol/gcdw/h')
        with self.pamodel:
            # change glucose uptake rate
            self.pamodel.change_reaction_bounds(rxn_id=self.substrate_uptake_id,
                                                lower_bound=-substrate, upper_bound=-substrate)
            # solve the model
            sol_pam = self.pamodel.optimize()

            if self.pamodel.solver.status == 'optimal' and self.pamodel.objective.value != 0:
                self.substrate_range += [substrate]
                self.fluxes.append(sol_pam.fluxes)
                # save the flux allocation coefficients
                flux_coeff = self.pamodel.flux_allocation_coefficients
                # we are only interested in the enzymes
                flux_coeff = flux_coeff[flux_coeff['constraint'] == 'enzyme'].set_index('enzyme_id')[
                    'coefficient'].to_frame().T
                flux_coeff['bin'] = index
                flux_coeff['substrate'] = substrate

                self.fac_df = pd.concat([self.fac_df, flux_coeff])

    def calculate_r_squared(self, bin, index, substrate):
        error = []
        # get reference datapoints
        # lower than upper bound
        ref_data_df = self.valid_data_df[self.valid_data_df[self.substrate_uptake_id] <= -bin[0]]
        #higher than lower bound
        ref_data_df = ref_data_df[ref_data_df[self.substrate_uptake_id] >= -bin[1]]

        # calculate error for different exchange rates
        for rxn in self.reactions_with_data + self.growth_rate:

            # only select the rows which are filled with data
            ref_data_rxn = ref_data_df.dropna(axis=0, subset=rxn)
            #if there are no reference data points, continue to the next reaction
            if len(ref_data_rxn) == 0: continue

            # calculate difference between simulations and validation data
            # get the linear relation of simulation (using the first and the last datapoint)
            line = linregress(x=[-b for b in bin[:2]], y=[self.fluxes[0][rxn], self.fluxes[-1][rxn]])
            ref_data_rxn = ref_data_rxn.assign(simulation=lambda x: line.intercept + line.slope * x[self.substrate_uptake_id])

            # simulation mean
            data_average = ref_data_rxn[rxn].mean()
            # error: squared difference
            ref_data_rxn = ref_data_rxn.assign(error=lambda x: (x[rxn] - x['simulation']) ** 2)

            # calculate R^2:
            residual_ss = sum(ref_data_rxn.error)
            total_ss = sum([(data - data_average) ** 2 for data in ref_data_rxn[rxn]])
            r_squared = 1 - residual_ss / total_ss
            error += [r_squared]
        self.error_df.loc[len(self.error_df)] = [index] + [substrate] + error
        print(self.error_df.to_markdown())

    def compare_sensitivities(self):
        #empty dataframes
        self.sensitive_enzymes = pd.DataFrame(columns=['bin', 'mean_sensitivity', 'enzyme_id'])
        self.bins_to_change = pd.DataFrame(columns=['bin', 'split', 'merge'])

        for i in self.bins:
            # get FAC data from bin
            fac_df_bin = self.fac_df[self.fac_df['bin'] == i]
            # if there is no solution, continue to the next bin
            if len(fac_df_bin) == 0:
                continue

            # 1. determine top n values based on average FAC
            fac_df_bin_t = fac_df_bin.drop(['bin', 'substrate'], axis=1).T

            fac_df_bin_t['fac_average'] = fac_df_bin_t.mean(axis=1)

            # convert the values to absolute to get the absolute highest values
            fac_df_bin_t['fac_average_absolute'] = [abs(avg) for avg in fac_df_bin_t['fac_average']]

            #sort the FAC and get the topn kcats to adjust
            fac_df_bin_t = fac_df_bin_t.sort_values(by='fac_average_absolute', ascending=False)
            fac_topn = fac_df_bin_t[:self.nmbr_kcats_to_change]
            # get ids and select those from the dataframe
            topn_ids = fac_topn.index.to_list()
            fac_df_topn = fac_df_bin[topn_ids + ['substrate']]
            fac_df_topn = fac_df_topn.sort_values(by='substrate')

            # 2. determine fraction of change as a measure of variability within the bin (fac_start/fac_end)
            for id in topn_ids:
                #check if the FACs are nonzero
                if fac_df_topn[id].iloc[0] != 0 and fac_df_topn[id].iloc[-1] != 0:
                    frac_fac_change = fac_df_bin[id].iloc[-1] / fac_df_topn[id].iloc[0]
                #if one is zero use the distance
                elif fac_df_topn[id].iloc[-1] == 0:
                    frac_fac_change = fac_df_topn[id].iloc[0]
                #else the fraction of change is 0
                else:
                    frac_fac_change = fac_df_topn[id].iloc[-1]
                # change can be end>start (<= 0.9) or end<start (>=1.1)
                if frac_fac_change <= 1 - self.change_threshold or frac_fac_change >= 1 + self.change_threshold:
                    self.bins_to_change.loc[len(self.bins_to_change)] = [i, True, False]

                # 3. Save sensitivity
                mean_sensitivity = fac_df_bin_t.loc[id]['fac_average']
                self.sensitive_enzymes.loc[len(self.sensitive_enzymes)] = [i, mean_sensitivity, id]

    def sample_duplicate_enzymes(self):
        # 4. Find enzymes with multiple entries. Pick one entry randomly for further processing
        duplicates = self.sensitive_enzymes[self.sensitive_enzymes.duplicated(['enzyme_id'])]
        self.sensitive_enzymes = self.sensitive_enzymes.drop_duplicates(subset=['enzyme_id'], keep=False)
        duplicate_enzymes = duplicates['enzyme_id'].drop_duplicates()
        for i, enz_id in duplicate_enzymes.items():
            # sample 1 row
            sampled_enz = duplicates[duplicates['enzyme_id'] == enz_id].sample(n=1)
            self.sensitive_enzymes = pd.concat([self.sensitive_enzymes, sampled_enz])

    def adjust_kcats_gradientdescent(self):
        for i, row in self.sensitive_enzymes.iterrows():
            enz_id = row['enzyme_id']
            fac = row['mean_sensitivity']
            rxn2kcat = self.pamodel.enzymes.get_by_id(enz_id).rxn2kcat
            new_rxn2kcat = {}
            for rxn, kcatdict in rxn2kcat.items():
                new_kcatdict = {}
                for direction, kcat in kcatdict.items():
                    new_kcat = kcat + kcat * fac * self.LEARNING_RATE
                    if new_kcat == 0.0 or np.isnan(new_kcat):
                        # print(new_kcat)
                        continue
                    else:
                        new_kcatdict[direction] = new_kcat
                if len(new_kcatdict) > 0:
                    new_rxn2kcat[rxn] = new_kcatdict
            if len(new_rxn2kcat) > 0:
                self.pamodel.change_kcat_value(enz_id, new_rxn2kcat)


#Setting the relative paths
BASE_DIR = os.path.split(os.getcwd())[0]
DATA_DIR = os.path.join(BASE_DIR, 'Data')
VALID_DATA_DIR = os.path.join(DATA_DIR, 'Ecoli_phenotypes')
valid_data_file = os.path.join(VALID_DATA_DIR, 'Ecoli_phenotypes_py_rev.xls')

valid_data_df = pd.read_excel(valid_data_file,sheet_name='Yields')

params = Parametrizer(valid_data_df, rxn_id)
params.run()