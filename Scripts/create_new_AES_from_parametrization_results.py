import pandas as pd
import os

from PAModelpy import CatalyticEvent

def search_index_in_parameter_file(df:pd.DataFrame, protein:str, reaction:str, direction:str):
    all_protein_rows = df[df.uniprot_id == protein]
    all_rxn_and_protein_rows = all_protein_rows[all_protein_rows.rxn_id == reaction]
    the_row = all_rxn_and_protein_rows[all_rxn_and_protein_rows.direction == direction]
    return the_row.index

RESULT_PARAMETRIZATION_FILE = os.path.join('Results', 'pam_parametrizer_diagnostics_ecolicore_UE0.7.xlsx')
NEW_AES_SUFFIX = 'ecolicore_240906'
SECTOR_PARAM_FILE = os.path.join('Data','proteinAllocationModel_iML1515_EnzymaticData_240730.xlsx')


if __name__ == '__main__':
    aes_parameter_file = pd.read_excel(SECTOR_PARAM_FILE, sheet_name='ActiveEnzymes')
    tps_parameter_file = pd.read_excel(SECTOR_PARAM_FILE, sheet_name='Translational')
    ues_parameter_file = pd.read_excel(SECTOR_PARAM_FILE, sheet_name='UnusedEnzyme')

    parametrization_results = pd.read_excel(RESULT_PARAMETRIZATION_FILE, sheet_name= 'Best_Individuals')

    #extract reaction id from catalytic event id
    parametrization_results['rxn_id'] = [CatalyticEvent._extract_reaction_id_from_catalytic_reaction_id(id) for id in parametrization_results['rxn_id']]
    #split all enzyme complex ids and make seperate rows from them
    parametrization_results['enzyme_id'] = parametrization_results.enzyme_id.str.split('_')
    parametrization_results = parametrization_results.explode('enzyme_id')

    for index, row in parametrization_results.iterrows():
        aes_index = search_index_in_parameter_file(aes_parameter_file, row.enzyme_id, row.rxn_id, row.direction)
        aes_parameter_file.loc[aes_index, 'kcat_values'] = row['kcat[s-1]']

    write_mode = 'w'
    kwargs = {}
    result_file = os.path.join('Results','parameter_files',f'proteinAllocationModel_EnzymaticData_{NEW_AES_SUFFIX}.xlsx')
    if os.path.isfile(result_file):
        write_mode = 'a'
        kwargs = {'if_sheet_exists':'replace'}

    with pd.ExcelWriter(os.path.join('Results','parameter_files',f'proteinAllocationModel_EnzymaticData_{NEW_AES_SUFFIX}.xlsx'),
                                     mode = write_mode, engine='openpyxl', **kwargs) as writer:
        aes_parameter_file.to_excel(writer, sheet_name = 'ActiveEnzymes', index=False)
        tps_parameter_file.to_excel(writer, sheet_name='Translational', index=False)
        ues_parameter_file.to_excel(writer, sheet_name='UnusedEnzyme', index=False)

