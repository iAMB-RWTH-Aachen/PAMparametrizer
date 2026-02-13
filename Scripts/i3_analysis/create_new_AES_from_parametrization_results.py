import pandas as pd
import os
import numpy as np

from PAModelpy import CatalyticEvent
from Modules.PAMparametrizer.utils.pamparametrizer_setup import save_sector_information_to_excel


RESULT_PARAMETRIZATION_FILE = os.path.join('Results', '2_parametrization', 'pam_parametrizer_diagnostics_iML1515_randomized.xlsx')
NEW_AES_SUFFIX = 'iMl1515_241209'
SECTOR_PARAM_FILE = os.path.join('Results','2_parametrization','proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx')


def search_index_in_parameter_file(df:pd.DataFrame, protein:str, reaction:str, direction:str):
    all_protein_rows = df[df.enzyme_id == protein]
    all_rxn_and_protein_rows = all_protein_rows[all_protein_rows.rxn_id == reaction]
    the_row = all_rxn_and_protein_rows[all_rxn_and_protein_rows.direction == direction]
    return the_row.index

def create_new_aes_parameter_file(old_param_file:str = SECTOR_PARAM_FILE,
                                  result_file_path:str = RESULT_PARAMETRIZATION_FILE,
                                  new_aes_suffix: str = NEW_AES_SUFFIX,
                                  default_enzyme_id_pattern: str = r'E[0-9][0-9]*|Enzyme_*'
                                  ) -> str:

    parameter_files = pd.read_excel(old_param_file, sheet_name=None)
    aes_parameter_file =parameter_files['ActiveEnzymes']

    parametrization_results = pd.read_excel(
        result_file_path, sheet_name='Best_Individuals').drop_duplicates(
        ['rxn_id', 'direction', 'enzyme_id'], keep='last')

    # extract reaction id from catalytic event id
    parametrization_results['rxn_id'] = [
        CatalyticEvent._extract_reaction_id_from_catalytic_reaction_id(id,
                                                                       default_enzyme_id_pattern = default_enzyme_id_pattern
                                                                       )
        for id in parametrization_results['rxn_id']
    ]

    for index, row in parametrization_results.iterrows():
        aes_index = search_index_in_parameter_file(aes_parameter_file, row.enzyme_id, row.rxn_id, row.direction)
        aes_parameter_file.loc[aes_index, 'kcat_values'] = row['kcat[s-1]']

    write_mode = 'w'
    kwargs = {}
    result_file = os.path.join('Results', '3_analysis', 'parameter_files',
                               f'proteinAllocationModel_EnzymaticData_{new_aes_suffix}.xlsx')
    if os.path.isfile(result_file):
        write_mode = 'a'
        kwargs = {'if_sheet_exists': 'replace'}

    with pd.ExcelWriter(result_file,
            mode=write_mode, engine='openpyxl', **kwargs) as writer:
        aes_parameter_file.to_excel(writer, sheet_name='ActiveEnzymes', index=False)
        for sheet, df in parameter_files.items():
            if sheet != 'ActiveEnzymes':
                df.to_excel(writer, sheet_name = sheet, index=False)

    return result_file

def change_enzyme_sector_in_excel(result_file_path: str,
                                  output_file_path: str,
                                  carbon_source: str,
                                  enzyme_sector: str = 'UnusedEnzymeSector') -> None:
    sectors_to_sheet = {'UnusedEnzymeSector':'UnusedEnzyme', 'TranslationalProteinSector': 'Translational'}
    parametrization_results = pd.read_excel(
        result_file_path, sheet_name='sector_parameters')
    sector_params = parametrization_results.loc[
        ((parametrization_results.substrate_uptake_id == carbon_source) &
         (parametrization_results.sector_id == enzyme_sector))
    ].rename({'substrate_uptake_id': 'lin_rxn_id'},
             axis=1)[['slope', 'intercept']].to_dict('records')[0]

    save_sector_information_to_excel(param_vs_lin_rxn = sector_params,
                                     lin_rxn_id=carbon_source,
                                     sector_id=sectors_to_sheet[enzyme_sector],
                                     pam_data_file = output_file_path,
                                     output_file_path = output_file_path
                                     )

def gaussian(x, mean, amplitude, standard_deviation):
    return amplitude * np.exp( - (x - mean)**2 / (2*standard_deviation ** 2))

def main_icgb21fr():
    new_ues_files = np.arange(1, 6, 1)
    for file_nmbr in range(1, 6):
        suffix = f'iCGB21FR_{file_nmbr}'
        result_file = os.path.join('Results', '2_parametrization', 'diagnostics',
                                   f'pam_parametrizer_diagnostics_{suffix}.xlsx')
        output_file_path = create_new_aes_parameter_file(
            old_param_file=os.path.join(
                'Results', '2_parametrization', 'proteinAllocationModel_iCGB21FR_EnzymaticData_multi.xlsx'
            ),
            result_file_path=result_file,
            new_aes_suffix=suffix,
            default_enzyme_id_pattern=r'Enzyme_cg[0-9]+'
        )

def main_ijn1463():
    for file_nmbr in range(1, 6):
        suffix = f'iJN1463_{file_nmbr}'
        result_file = os.path.join('Results', '2_parametrization', 'diagnostics',
                                   f'pam_parametrizer_diagnostics_{suffix}.xlsx')
        create_new_aes_parameter_file(
            old_param_file=os.path.join(
                'Results', '2_parametrization', 'proteinAllocationModel_iJN1463_EnzymaticData_multi.xlsx'
            ),
            result_file_path=result_file,
            new_aes_suffix=suffix,
            default_enzyme_id_pattern=r'Enzyme_*|Enzyme_PP_[0-9]+'
        )


def main_iml1515():
    new_ues_files = [1, 2, 4, 5, 6]
    for file_nmbr in range(1, 11):
        # suffix = f'iML1515_{file_nmbr}'
        suffix = 'iML1515_csources'
        result_file = os.path.join('Results', '2_parametrization', 'diagnostics',
                                   f'pam_parametrizer_diagnostics_{suffix}.xlsx')
        output_file_path = create_new_aes_parameter_file(
            old_param_file=os.path.join(
                'Results', '2_parametrization', 'proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx'
            ),
            result_file_path=result_file,
            new_aes_suffix=suffix,
        )

        if file_nmbr in new_ues_files:
            change_enzyme_sector_in_excel(result_file_path=result_file,
                                          output_file_path=output_file_path,
                                          carbon_source='EX_glc__D_e')

            change_enzyme_sector_in_excel(result_file_path=result_file,
                                          output_file_path=output_file_path,
                                          enzyme_sector='TranslationalProteinSector',
                                          carbon_source='EX_glc__D_e')

if __name__ == '__main__':
    # main_ijn1463()
    main_iml1515()

