import os
import pandas as pd
import numpy as np

from Modules.utils.pam_generation import setup_pputida_pam
from Modules.utils.pam_generation import create_pamodel_from_diagnostics_file
from Modules.utils.pamparametrizer_analysis import get_results_from_simulations

NUM_MODELS = 5

PAM_KCAT_FILES_IJN = [os.path.join('Results', '2_parametrization', 'diagnostics',
                               f'pam_parametrizer_diagnostics_iJN1463_{file_nmbr}.xlsx') for file_nmbr in
                  range(1, NUM_MODELS + 1)]

PPUTIDA_PHENOTYPE_FILE_PATH = os.path.join('Data', 'Pputida_phenotypes', 'pputida_phenotypes.xlsx')
RXNS_TO_VALIDATE = {'Peripheral': ['GNK'],
                    'EMP': ['PGI', 'ENO'],
                    'PPP': ['TKT2'],
                    'ED': ['EPA'],
                    'TCA': ['PDH', 'SUCDi']
                    }
rxns_to_save = []
for rxns in RXNS_TO_VALIDATE.values(): rxns_to_save+=rxns

if __name__ == '__main__':
    mfa_data = pd.read_excel(PPUTIDA_PHENOTYPE_FILE_PATH,
                             sheet_name='fluxomics_glucose').iloc[0] # first row has most relevant data
    pam = setup_pputida_pam(sensitivity=False)

    fluxes = {'GotEnzymes': get_results_from_simulations(pam,
                                                         substrate_rates = [[mfa_data['EX_glc__D_e']]],
                                                         fluxes_to_save=rxns_to_save,
                                                         transl_sector_config = False
                                               )['fluxes']
              }

    for i, alternative_file in enumerate(PAM_KCAT_FILES_IJN):
        alt_pam = create_pamodel_from_diagnostics_file(
            alternative_file,
            pam.copy(copy_with_pickle=True,
                     other_enzyme_id_pattern = r'E[0-9][0-9]*|Enzyme*|PP_*')
        )
        fluxes[f'Alternative {i}'] = get_results_from_simulations(alt_pam,
                                                                  substrate_rates = [[mfa_data['EX_glc__D_e']]],
                                                                  fluxes_to_save=rxns_to_save,
                                                                  transl_sector_config=False
                                               )['fluxes']
flux_df = pd.DataFrame()
for model, df in fluxes.items():
    df['model']=model
    flux_df = pd.concat([flux_df, df])

flux_df.to_excel('Results/3_analysis/mfa_experiment_pputida.xlsx')
print(flux_df.to_latex())



