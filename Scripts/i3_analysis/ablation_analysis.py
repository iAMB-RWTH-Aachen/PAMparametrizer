import pandas as pd
import numpy as np
from typing import List
import os

from cobra.io import read_sbml_model
from PAModelpy.utils import set_up_pam

from Modules.PAMparametrizer.utils.pam_generation import create_pamodel_from_diagnostics_file



PARAM_FILE_OLD = os.path.join('Results', '1_preprocessing','proteinAllocationModel_iML1515_EnzymaticData_250912.xlsx')
PARAM_FILE_PREPROC = os.path.join('Results','2_parametrization','proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx')
EXCHANGE_FLUX_FILE_PATH = os.path.join('Data', 'Ecoli_phenotypes', 'Ecoli_phenotypes_py.xls')
MFA_DATA_FILE_PATH = os.path.join('Data', 'Ecoli_phenotypes', 'fluxomics_datasets_gerosa.xlsx')

MODEL_FILE = os.path.join('Models', 'iML1515.xml')
SUBSTRATE_ID = 'EX_glc__D_e'
GEROSA_GLC_UPTAKE = -9.654

RESULT_FLUX_PATH = os.path.join('Results', '3_analysis', 'iML1515_alternative_models_predictions_diff_growth_rates.xlsx')

def save_flux_predictions(model_files: List[str],
                         substrate_uptake_rates: List[float],
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
        for substrate_uptake_rate in substrate_uptake_rates:

            model.reactions.get_by_id(SUBSTRATE_ID).lower_bound = substrate_uptake_rate
            model.optimize()

            all_reactions = [rxn for rxn in model.reactions if 'CE_' not in rxn.id]#CE to ignore the catalytic events
            fluxes = (pd.DataFrame.from_dict({rxn.id: rxn.flux for rxn in all_reactions},
                                  columns =['flux'], orient='index',)
                      .assign(substrate_uptake_rate = substrate_uptake_rate)
                      )

        kwargs = {'mode': 'a', 'if_sheet_exists': 'replace'} \
            if os.path.exists(simulated_fluxes_result_file) else {'mode': 'w'}
        with pd.ExcelWriter(simulated_fluxes_result_file,
                            engine='openpyxl', **kwargs) as writer:
            fluxes.to_excel(writer, sheet_name=label)

def perform_and_save_simulations(num_alternative_models: int,
                                 substrate_uptake_rates: List[float],
                                 simulated_fluxes_result_file: str = RESULT_FLUX_PATH,
                                 ):
    NUM_ALT_MODELS = 10
    diagnostic_files = list()
    param_files = list()
    all_model_labels = ['iML1515', 'GotEnzymes', 'After preprocessing'] \
                       + [f'alternative {i}' for i in range(1, NUM_ALT_MODELS + 1)]
    for file_nmbr in range(1, NUM_ALT_MODELS + 1):
        diagnostic_files += [os.path.join('Results', '2_parametrization', 'diagnostics',
                                          f'pam_parametrizer_diagnostics_{file_nmbr}.xlsx')]
        param_files += [os.path.join('Results', '3_analysis', 'parameter_files',
                                     f'proteinAllocationModel_EnzymaticData_iML1515_{file_nmbr}.xlsx')]

    save_flux_predictions(model_files=[MODEL_FILE, PARAM_FILE_OLD,
                                      PARAM_FILE_PREPROC],
                          substrate_uptake_rates=substrate_uptake_rates,
                          diagnostic_files=diagnostic_files,
                          label_names=all_model_labels,
                          simulated_fluxes_result_file=simulated_fluxes_result_file,
                          )

if __name__ == '__main__':
    mfa_data = pd.read_excel(MFA_DATA_FILE_PATH)
    mfa_data_glc = mfa_data[mfa_data.condition == 'Glucose'][['reaction', 'measured']].set_index('reaction').squeeze()

    exchange_fluxes = pd.read_excel(EXCHANGE_FLUX_FILE_PATH, sheet_name='Yields')
    exchange_fluxes = exchange_fluxes[exchange_fluxes.Strain == 'MG1655']

    print(mfa_data_glc)
    print(exchange_fluxes)