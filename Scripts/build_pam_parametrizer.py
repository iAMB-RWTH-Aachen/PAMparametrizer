import os
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")
from PAModelpy.configuration import Config

from Modules.PAM_parametrizer import ValidationData, HyperParameters, ParametrizationResults
from Modules.PAM_parametrizer import PAMParametrizer
from Scripts.pam_generation import setup_toy_pam, setup_ecolicore_pam

MAX_SUBSTRATE_UPTAKE_RATE_toy = 0.1
MIN_SUBSTRATE_UPTAKE_RATE_toy = 0.001
MAX_SUBSTRATE_UPTAKE_RATE_ecoli = -0.1
MIN_SUBSTRATE_UPTAKE_RATE_ecoli = -10



def set_up_validation_data():
    config = Config()
    config.reset()
    RXNS_TO_VALIDATE = [config.ACETATE_EXCRETION_RXNID, config.CO2_EXHANGE_RXNID, config.OXYGEN_UPTAKE_RXNID,
                        'BIOMASS_Ecoli_core_w_GAM']
    DATA_DIR = os.path.join(os.getcwd(), 'Data')
    VALID_DATA_PATH = os.path.join(DATA_DIR, 'Ecoli_phenotypes', 'Ecoli_phenotypes_py_rev.xls')

    valid_data_df = pd.read_excel(VALID_DATA_PATH,sheet_name='Yields')
    valid_data_df = valid_data_df.rename(columns = {config.GLUCOSE_EXCHANGE_RXNID: config.GLUCOSE_EXCHANGE_RXNID+'_ub',
                                                    config.BIOMASS_REACTION: 'BIOMASS_Ecoli_core_w_GAM'})
    validation_data = ValidationData(valid_data_df)
    validation_data._reactions_to_plot = RXNS_TO_VALIDATE
    validation_data._reactions_to_validate = RXNS_TO_VALIDATE
    return validation_data


def set_up_validation_data_toy_model():
    DATA_DIR = os.path.join(os.getcwd(), 'Scripts', 'Testing', 'Data')
    RESULT_DF_FILE = os.path.join(DATA_DIR, 'toy_model_simulations_ga.csv')
    valid_data_df = pd.read_csv(RESULT_DF_FILE)

    validation_data = ValidationData(valid_data_df)
    validation_data._reactions_to_plot = ['R1', 'R7', 'R8', 'R9']
    validation_data._reactions_to_validate = ['R1', 'R7', 'R8', 'R9']
    return validation_data

def set_up_hyperparameter():
    hyperparams = HyperParameters
    hyperparams.threshold_iteration = 10
    hyperparams.threshold_error = 0.95
    hyperparams.number_of_kcats_to_mutate = 4
    hyperparams.genetic_algorithm_hyperparams['number_generations'] = 6
    # hyperparams.genetic_algorithm_filename_base = 'genetic_algorithm_run_toy_'
    hyperparams.genetic_algorithm_hyperparams['print_progress'] = False
    return hyperparams

def set_up_toy_model(kcat_fwd:list = [1, 0.5, 1, 0.5, 0.45, 1.5]):
    toy_pam = setup_toy_pam(kcat_fwd=kcat_fwd)
    return toy_pam

def run_simulations(pamodel, substrate_rates:list,
                    reactions_to_eval:list, substrate_uptake_reaction:str):
    result_df = pd.DataFrame(columns= [substrate_uptake_reaction+'_ub']+reactions_to_eval)

    for substrate in substrate_rates:
        if substrate >0:
            pamodel.change_reaction_bounds(rxn_id=substrate_uptake_reaction,
                                       lower_bound=substrate, upper_bound=0)
        else:
            pamodel.change_reaction_bounds(rxn_id=substrate_uptake_reaction,
                                           lower_bound=0, upper_bound=substrate)
        print('Running simulations with ', substrate, 'mmol/g_cdw/h of substrate going into the system')
        pamodel.optimize()
        if pamodel.solver.status == 'optimal' and pamodel.objective.value>0:
            results_row = []
            for rxn in reactions_to_eval:
                results_row += [pamodel.reactions.get_by_id(rxn).flux]

            result_df.loc[len(result_df)] = [abs(substrate)] + results_row
    return result_df

def set_up_pamparametrizer(pamodel,
                           min_substrate_uptake_rate: float,
                           max_substrate_uptake_rate: float,
                           file_extension: str = ''):
    if pamodel.id == 'toy_model':
        validation_data = set_up_validation_data_toy_model()
    else:
        validation_data = set_up_validation_data()
    hyperparameters = set_up_hyperparameter()
    hyperparameters.genetic_algorithm_filename_base =  'genetic_algorithm_'+ file_extension

    return PAMParametrizer(pamodel=pamodel,
                     validation_data=validation_data,
                     hyperparameters=hyperparameters,
                     substrate_uptake_id='R1',
                     max_substrate_uptake_rate=max_substrate_uptake_rate,
                     min_substrate_uptake_rate=min_substrate_uptake_rate)

def run_pam_parametrizer_toymodel(binned:str = 'False',
                                  filename_extension: str =''):
    toy_pam = set_up_toy_model()
    pam_parametrizer = set_up_pamparametrizer(pamodel = toy_pam,
                                              min_substrate_uptake_rate=MIN_SUBSTRATE_UPTAKE_RATE_toy,
                                              max_substrate_uptake_rate=MAX_SUBSTRATE_UPTAKE_RATE_toy,
                                              file_extension= 'toy_'+filename_extension)
    pam_parametrizer.run(remove_subruns=True, binned = binned)
    return pam_parametrizer

def run_pam_parametrizer_ecolicore(binned:str = 'False',
                                  filename_extension: str =''):
    ecolicore_pam = setup_ecolicore_pam()
    pam_parametrizer = set_up_pamparametrizer(pamodel = ecolicore_pam,
                                              min_substrate_uptake_rate=MIN_SUBSTRATE_UPTAKE_RATE_ecoli,
                                              max_substrate_uptake_rate=MAX_SUBSTRATE_UPTAKE_RATE_ecoli,
                                              file_extension= 'ecolicore_'+ filename_extension)
    pam_parametrizer.run(remove_subruns=True, binned = binned)
    return pam_parametrizer

