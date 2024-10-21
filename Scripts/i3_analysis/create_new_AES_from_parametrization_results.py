import matplotlib.pyplot as plt
from matplotlib.colors import to_hex

import pandas as pd
import os
import numpy as np
from PAModelpy import CatalyticEvent


RESULT_PARAMETRIZATION_FILE = os.path.join('Results', '2_parametrization', 'pam_parametrizer_diagnostics_iML1515_randomized.xlsx')
NEW_AES_SUFFIX = 'iMl1515_241009'
SECTOR_PARAM_FILE = os.path.join('Data','proteinAllocationModel_iML1515_EnzymaticData_240730_multi.xlsx')
RESULT_HISTOGRAM_PATH = os.path.join('Results', '3_analysis', 'aes_histogram_iML1515_241012.png')

def search_index_in_parameter_file(df:pd.DataFrame, protein:str, reaction:str, direction:str):
    all_protein_rows = df[df.uniprot_id == protein]
    all_rxn_and_protein_rows = all_protein_rows[all_protein_rows.rxn_id == reaction]
    the_row = all_rxn_and_protein_rows[all_rxn_and_protein_rows.direction == direction]
    return the_row.index

def create_new_aes_parameter_file(old_param_file:str = SECTOR_PARAM_FILE,
                                  result_file_path:str = RESULT_PARAMETRIZATION_FILE,
                                  new_aes_suffix: str = NEW_AES_SUFFIX):
    aes_parameter_file = pd.read_excel(old_param_file, sheet_name='ActiveEnzymes')
    tps_parameter_file = pd.read_excel(old_param_file, sheet_name='Translational')
    ues_parameter_file = pd.read_excel(old_param_file, sheet_name='UnusedEnzyme')

    parametrization_results = pd.read_excel(result_file_path, sheet_name='Best_Individuals')

    # extract reaction id from catalytic event id
    parametrization_results['rxn_id'] = [CatalyticEvent._extract_reaction_id_from_catalytic_reaction_id(id) for id in
                                         parametrization_results['rxn_id']]
    # split all enzyme complex ids and make seperate rows from them
    parametrization_results['enzyme_id'] = parametrization_results.enzyme_id.str.split('_')
    parametrization_results = parametrization_results.explode('enzyme_id')

    for index, row in parametrization_results.iterrows():
        aes_index = search_index_in_parameter_file(aes_parameter_file, row.enzyme_id, row.rxn_id, row.direction)
        aes_parameter_file.loc[aes_index, 'kcat_values'] = row['kcat[s-1]']

    write_mode = 'w'
    kwargs = {}
    result_file = os.path.join('Results', '1_preprocessing',
                               f'proteinAllocationModel_EnzymaticData_{new_aes_suffix}.xlsx')
    if os.path.isfile(result_file):
        write_mode = 'a'
        kwargs = {'if_sheet_exists': 'replace'}

    with pd.ExcelWriter(
            os.path.join('Results', '1_preprocessing', f'proteinAllocationModel_EnzymaticData_{new_aes_suffix}.xlsx'),
            mode=write_mode, engine='openpyxl', **kwargs) as writer:
        aes_parameter_file.to_excel(writer, sheet_name='ActiveEnzymes', index=False)
        tps_parameter_file.to_excel(writer, sheet_name='Translational', index=False)
        ues_parameter_file.to_excel(writer, sheet_name='UnusedEnzyme', index=False)

def gaussian(x, mean, amplitude, standard_deviation):
    return amplitude * np.exp( - (x - mean)**2 / (2*standard_deviation ** 2))

def create_kcat_histogram_old_vs_new(data_file_paths: list[pd.DataFrame],
                                     label_names:list[str],
                                     result_fig_file: str = RESULT_HISTOGRAM_PATH):
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
                                              stacked=True, fill=False, label= label, color=color)
        i+=1
        peak_bin_index = np.argmax(bin_heights)
        peak_bin_value = (bin_borders[peak_bin_index]+bin_borders[peak_bin_index+1])/2

        print(f'\tMost frequent kcat:\t{peak_bin_value}')

    ax.vlines([13.7], 0, 1e4, linestyles='dotted')

    plt.yscale('log')
    plt.ylabel('Frequency')
    plt.xscale('log')
    plt.xlabel('Kcat value [s-1]')

    plt.legend()
    plt.tight_layout()
    plt.savefig(result_fig_file)

if __name__ == '__main__':
    other_files = [os.path.join('Results', '1_preprocessing',
                               f'proteinAllocationModel_EnzymaticData_iMl1515_241009.xlsx')]
    for file_nmbr in [1,2]:
        suffix = f'iML1515_{file_nmbr}'
        # result_file = os.path.join('Results', f'pam_parametrizer_diagnostics_{file_nmbr}.xlsx')
        # create_new_aes_parameter_file(result_file_path= result_file,
        #                               new_aes_suffix= suffix)
        other_files += [os.path.join('Results', '1_preprocessing',
                               f'proteinAllocationModel_EnzymaticData_{suffix}.xlsx')]

    create_kcat_histogram_old_vs_new([os.path.join('Data','proteinAllocationModel_iML1515_EnzymaticData_240730.xlsx'),
                                      SECTOR_PARAM_FILE] + other_files,
                                     label_names = ['GotEnzymes', 'After preprocessing']\
                                                   + [f'alternative {i}' for i in [1,2,3]])

