import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import numpy as np
import pandas as pd
import os

from Scripts.pam_generation_uniprot_id import set_up_ecoli_pam

SECTOR_PARAM_FILE = os.path.join('Data','proteinAllocationModel_iML1515_EnzymaticData_240730_multi.xlsx')
RESULT_KCAT_HISTOGRAM_PATH = os.path.join('Results', '3_analysis', 'aes_histogram_iML1515_241012.png')
RESULT_FLUX_HISTOGRAM_PATH = os.path.join('Results', '3_analysis', 'flux_histogram_iML1515_241012.png')

def gaussian(x, mean, amplitude, standard_deviation):
    return amplitude * np.exp( - (x - mean)**2 / (2*standard_deviation ** 2))

def calculate_distribution_statistics(bin_heigths: list[float],
                                      bin_borders:list[float]) -> None:
    peak_bin_index = np.argmax(bin_heigths)
    peak_bin_value = (bin_borders[peak_bin_index] + bin_borders[peak_bin_index + 1]) / 2
    area_under_curve = sum(
        [(abs(bin_borders[i]) - abs(bin_borders[i + 1])) * bin_heigths[i] for i in range(len(bin_heigths))])

    print(
        f'\tMost frequent flux:\t\t{peak_bin_value} mmol/gCDW/h\n\tArea under the curve:\t\t{area_under_curve} mmol/gCDW/h\n')

def create_kcat_histogram_old_vs_new(data_file_paths: list[pd.DataFrame],
                                     label_names:list[str],
                                     result_fig_file: str = RESULT_KCAT_HISTOGRAM_PATH,
                                     cumulative: bool = False):
    fig, ax = plt.subplots()
    n_bins = 50
    i = 0
    cmap = plt.get_cmap("viridis")

    for label, data_file_path in zip(label_names, data_file_paths):
        aes_parameter_df = pd.read_excel(data_file_path, sheet_name='ActiveEnzymes')
        kcat_values = aes_parameter_df.kcat_values.dropna()
        print('------------------------------------------------------------------------')
        print(f'The kcat set from {label} has:\n \tMedian:\t\t\t{np.median(kcat_values)} \n \tMean:\t\t\t{np.mean(kcat_values)}')


        hist, bins = np.histogram(kcat_values, bins=n_bins)
        logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))

        color = to_hex(cmap(i / (len(data_file_paths))))
        bin_heights, bin_borders, _ = ax.hist(kcat_values, bins = logbins, histtype='step',
                                              stacked=True, fill=False, label= label, color=color, cumulative=cumulative)
        i+=1
        calculate_distribution_statistics(bin_heights, bin_borders)

    ax.vlines([13.7], 0, 1e4, linestyles='dotted')

    plt.yscale('log')
    plt.ylabel('Frequency')
    plt.xscale('log')
    plt.xlabel('Kcat value [s-1]')

    plt.legend()
    plt.tight_layout()
    plt.savefig(result_fig_file)

def create_flux_histogram_old_vs_new(data_file_paths: list[pd.DataFrame],
                                     label_names:list[str],
                                     result_fig_file: str = RESULT_FLUX_HISTOGRAM_PATH,
                                     cumulative:bool=False):
    fig, ax = plt.subplots()
    n_bins = 50
    i = 0
    cmap = plt.get_cmap("viridis")

    for label, data_file_path in zip(label_names, data_file_paths):
        print('------------------------------------------------------------------------')
        print('Analyzing flux distribution of model', label, '\n')
        model = set_up_ecoli_pam(pam_info_file=data_file_path, sensitivity=False)
        model.change_reaction_bounds('EX_glc__D_e', -11,0)
        solution = model.optimize()

        fluxes = [abs(flux) for flux in solution.fluxes.values if flux!=0]
        print(
            f'Model {label} has:\n \ta growth rate of:\t\t{model.objective.value} h-1 with 11 mmol_glc/gCDW/h \n '
            f'\tMedian fluxes:\t\t\t{np.median(fluxes)} mmol/CgDW/h\n\tMean fluxes:\t\t\t{np.mean(fluxes)} mmol/gCDW/h')

        hist, bins = np.histogram(fluxes, bins=n_bins)
        logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))

        color = to_hex(cmap(i / (len(data_file_paths))))
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
    other_files = [os.path.join('Results', '3_analysis', 'parameter_files',
                               'proteinAllocationModel_EnzymaticData_iML1515_241009.xlsx')]

    for file_nmbr in [1,2]:
        suffix = f'iML1515_{file_nmbr}'

        other_files += [os.path.join('Results', '3_analysis', 'parameter_files',
                               f'proteinAllocationModel_EnzymaticData_{suffix}.xlsx')]

    create_flux_histogram_old_vs_new([os.path.join('Data','proteinAllocationModel_iML1515_EnzymaticData_240730.xlsx'),
                                      SECTOR_PARAM_FILE] + other_files,
                                     label_names = ['GotEnzymes', 'After preprocessing']\
                                                   + [f'alternative {i}' for i in [1,2,3]])
    create_kcat_histogram_old_vs_new([os.path.join('Data', 'proteinAllocationModel_iML1515_EnzymaticData_240730.xlsx'),
                                      SECTOR_PARAM_FILE] + other_files,
                                     label_names=['GotEnzymes', 'After preprocessing'] \
                                                 + [f'alternative {i}' for i in [1, 2, 3]])
