import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import seaborn as sns
import numpy as np
import pandas as pd
import os

from PAModelpy.utils import set_up_pam
from Figures.Figure1_iml1515_kcat_analysis import calculate_distribution_statistics


MODEL_FILE = os.path.join('Models', 'iML1515.xml')
SUBSTRATE_ID = 'EX_glc__D_e'

def create_flux_histogram_old_vs_new(data_file_paths: list[pd.DataFrame],
                                     label_names:list[str],
                                     result_fig_file: str,
                                     other_colors = {'GotEnzymes': 'grey', 'After preprocessing': 'black'},
                                     cumulative:bool=False):
    fig, ax = plt.subplots()
    n_bins = 50
    i = 0
    cmap = plt.get_cmap("coolwarm")

    for label, data_file_path in zip(label_names, data_file_paths):
        print('------------------------------------------------------------------------')
        print('Analyzing flux distribution of model', label, '\n')
        model = set_up_pam(pam_info_file=data_file_path,
                           model=MODEL_FILE,
                           sensitivity=False)
        model.change_reaction_bounds(SUBSTRATE_ID, -11,0)
        model.optimize()

        fluxes = [abs(rxn.flux) for rxn in model.reactions if rxn.flux!=0 and 'CE_' not in rxn.id] #CE to ignore the catalytic events
        print(
            f'Model {label} has:\n \ta growth rate of:\t\t{model.objective.value} h-1 with 11 mmol_glc/gCDW/h \n '
            f'\tMedian fluxes:\t\t\t{np.median(fluxes)} mmol/gCDW/h\n\tMean fluxes:\t\t\t{np.mean(fluxes)} mmol/gCDW/h')

        hist, bins = np.histogram(fluxes, bins=n_bins)
        logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))

        if label in other_colors.keys():
            color = other_colors[label]
        else:
            i += 1
            color = to_hex(cmap(i / (len(data_file_paths) - len(other_colors))))
        kwargs = {'cumulative':cumulative}
        if not cumulative:
            kwargs = {**kwargs, 'fill': True,  'alpha':0.5}

        bin_heights, bin_borders, _ = ax.hist(fluxes, bins=logbins, histtype='step',
                                              stacked = True, label = label, color = color,
                                              **kwargs)

        i += 1
        calculate_distribution_statistics(bin_heights, bin_borders)


    # plt.yscale('log')
    plt.ylabel('Frequency')
    plt.xscale('log')
    plt.xlabel('Flux [mmol/gCDW/h]')

    plt.legend()
    plt.tight_layout()
    plt.savefig(result_fig_file)

if __name__ == '__main__':
    NUM_ALT_MODELS = 8
    PARAM_FILE_ORI = os.path.join('Results', '1_preprocessing',
                                  'proteinAllocationModel_iML1515_EnzymaticData_250225.xlsx')
    PARAM_FILE_PREPROC = os.path.join('Results', '2_parametrization',
                                      'proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx')

    other_files = list()
    for file_nmbr in range(1, NUM_ALT_MODELS + 1):
        suffix = f'iML1515_{file_nmbr}'

        other_files += [os.path.join('Results', '3_analysis', 'parameter_files',
                                     f'proteinAllocationModel_EnzymaticData_{suffix}.xlsx')]

    create_flux_histogram_old_vs_new([PARAM_FILE_ORI,
                                      PARAM_FILE_PREPROC] + other_files,
                                     label_names=['GotEnzymes', 'After preprocessing'] \
                                                 + [f'alternative {i}' for i in range(1, NUM_ALT_MODELS + 1)],
                                     cumulative=True)