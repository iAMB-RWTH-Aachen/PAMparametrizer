from cobra.io.sbml import read_sbml_model
import pandas as pd
import os
from typing import Union, List
from PAModelpy import PAModel
from PAModelpy.utils import set_up_pam
from cobra.flux_analysis import flux_variability_analysis
from Modules.PAMparametrizer.utils.pam_generation import create_pamodel_from_diagnostics_file


PARAM_FILE_OLD = os.path.join('Results', '1_preprocessing','proteinAllocationModel_iML1515_EnzymaticData_250912.xlsx')
PARAM_FILE_PREPROC = os.path.join('Results','2_parametrization','proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx')
MODEL_FILE = os.path.join('Models', 'iML1515.xml')
SUBSTRATE_ID = 'EX_glc__D_e'
GEROSA_GLC_UPTAKE = -9.654

RESULT_FLUX_PATH = os.path.join('Results', '3_analysis', 'iML1515_alternative_models_preditions.xlsx')


def get_all_kcat_values(data_file_paths: list[pd.DataFrame],
                                     label_names:list[str]):


    all_kcat_values = pd.DataFrame()
    for label, data_file_path in zip(label_names,data_file_paths):
        aes_parameter_df = pd.read_excel(data_file_path, sheet_name='ActiveEnzymes')[['rxn_id', 'enzyme_id', 'direction', 'kcat_values']]
        all_kcat_values = pd.merge(all_kcat_values, aes_parameter_df, on = ['rxn_id', 'enzyme_id', 'direction'],
                                    how ='outer', suffixes=['', label])
    return all_kcat_values

def get_flux_predictions_and_fva(model_files: List[str],
                                 diagnostic_files: List[str],
                                     label_names:List[str],
                                     simulated_fluxes_result_file: str = RESULT_FLUX_PATH,
                                     ):
    pam = set_up_pam(pam_info_file = PARAM_FILE_PREPROC, model = MODEL_FILE, sensitivity = False)
    models = [set_up_pam(pam_info_file = file, model = MODEL_FILE, sensitivity = False)
              if 'xml' not in file else read_sbml_model(file)
              for file in model_files
    ]
    alt_pams = [create_pamodel_from_diagnostics_file(diag_file, pam.copy_with_pickle()) for diag_file in diagnostic_files]

    for label, model in zip(label_names, models+alt_pams):
        print('------------------------------------------------------------------------')
        print('Analyzing flux distribution of model', label, '\n')

        model.reactions.get_by_id(SUBSTRATE_ID).lower_bound = GEROSA_GLC_UPTAKE
        model.optimize()

        all_reactions = [rxn for rxn in model.reactions if 'CE_' not in rxn.id]#CE to ignore the catalytic events
        fluxes = pd.DataFrame({rxn: rxn.flux for rxn in all_reactions},
                              columns =['rxn_id', 'flux']).set_index('rxn_id')
        fva_results = flux_variability_analysis(model,all_reactions,)
        all_results = pd.merge(fluxes, fva_results, left_index=True, right_index=True, how='outer')
        with pd.ExcelWriter(simulated_fluxes_result_file, engine='openpyxl',
                            mode='a', if_sheet_exists='replace') as writer:
            all_results.to_excel(writer, sheet_name=label, index=False)

def determine_flux_coefficient_of_variation():
    fva_results = pd.read_excel(RESULT_FLUX_PATH, sheet_name=None)
    parsed_fva_results ={}
    for model_id, df in fva_results.items():
        df['flux_cv'] = df.apply(lambda row: (row['max'] - row['min'])/row['flux'], axis=1)
        parsed_fva_results[model_id] = df
    return parsed_fva_results

def determine_kcat_coefficient_of_variation(kcat_df:pd.DataFrame):
    kcat_df['kcat_cv'] = kcat_df.apply(lambda row: row.std()/row.mean(), axis=1)
    return kcat_df

if __name__ == '__main__':
    NUM_ALT_MODELS = 10
    diagnostic_files = list()
    for file_nmbr in range(1, NUM_ALT_MODELS + 1):
        diagnostic_files += [os.path.join('Results', '2_parametrization', 'diagnostics',
                                     f'pam_parametrizer_diagnostics_{file_nmbr}.xlsx')]
    get_flux_predictions_and_fva(model_files=[MODEL_FILE, PARAM_FILE_OLD,
                                      PARAM_FILE_PREPROC], diagnostic_files=diagnostic_files,
                                     label_names=['iML1515','GotEnzymes', 'After preprocessing'] \
                                                 + [f'alternative {i}' for i in range(1, NUM_ALT_MODELS + 1)])



