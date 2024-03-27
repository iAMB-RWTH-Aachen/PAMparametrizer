import os
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")



from PAModelpy.configuration import Config

from Modules.PAM_parametrizer import ValidationData, HyperParameters, ParametrizationResults
from Modules.PAM_parametrizer import PAMParametrizer
from Scripts.pam_generation import setup_ecolicore_pam

# import sys
# sys.stdout = open('output.txt','wt')

MAX_SUBSTRATE_UPTAKE_RATE = -0.1
MIN_SUBSTRATE_UPTAKE_RATE = -10
config = Config()
config.reset()
RXNS_TO_VALIDATE = [config.ACETATE_EXCRETION_RXNID, config.CO2_EXHANGE_RXNID, config.OXYGEN_UPTAKE_RXNID, 'BIOMASS_Ecoli_core_w_GAM']

def set_up_validation_data():
    DATA_DIR = os.path.join(os.getcwd(), 'Data')
    VALID_DATA_PATH = os.path.join(DATA_DIR, 'Ecoli_phenotypes', 'Ecoli_phenotypes_py_rev.xls')

    valid_data_df = pd.read_excel(VALID_DATA_PATH,sheet_name='Yields')
    valid_data_df = valid_data_df.rename(columns = {config.GLUCOSE_EXCHANGE_RXNID: config.GLUCOSE_EXCHANGE_RXNID+'_ub',
                                                    config.BIOMASS_REACTION: 'BIOMASS_Ecoli_core_w_GAM'})
    validation_data = ValidationData(valid_data_df)
    validation_data._reactions_to_plot = RXNS_TO_VALIDATE
    validation_data._reactions_to_validate = [RXNS_TO_VALIDATE[0], RXNS_TO_VALIDATE[-1]]
    return validation_data

def set_up_hyperparameter():
    hyperparams = HyperParameters
    hyperparams.threshold_iteration = 10
    hyperparams.number_of_kcats_to_mutate = 4
    hyperparams.filename_extension = 'ecolicore_before'
    hyperparams.genetic_algorithm_hyperparams['number_generations'] = 6
    hyperparams.genetic_algorithm_filename_base = 'genetic_algorithm_run_ecolicore_'
    hyperparams.genetic_algorithm_hyperparams['print_progress'] = True
    return hyperparams

def run_simulations(pamodel, substrate_rates, rxn_to_validate = RXNS_TO_VALIDATE):
    result_df = pd.DataFrame(columns= rxn_to_validate)

    for substrate in substrate_rates:
        pamodel.change_reaction_bounds(rxn_id=config.GLUCOSE_EXCHANGE_RXNID,
                                       lower_bound=substrate, upper_bound=0)
        print('Running simulations with ', substrate, 'mmol/g_cdw/h of substrate going into the system')
        pamodel.optimize()
        if pamodel.solver.status == 'optimal' and pamodel.objective.value>0:
            results_row = []
            for rxn in rxn_to_validate:
                results_row += [pamodel.reactions.get_by_id(rxn).flux]

            result_df.loc[len(result_df)] = [substrate] + results_row
    return result_df

def set_up_pamparametrizer(min_substrate_uptake_rate:float, max_substrate_uptake_rate: float):
    ecolicore_pam = setup_ecolicore_pam()
    validation_data = set_up_validation_data()
    hyperparameters = set_up_hyperparameter()

    return PAMParametrizer(pamodel=ecolicore_pam,
                     validation_data=validation_data,
                     hyperparameters=hyperparameters,
                     substrate_uptake_id=config.GLUCOSE_EXCHANGE_RXNID,
                     max_substrate_uptake_rate=max_substrate_uptake_rate,
                     min_substrate_uptake_rate=min_substrate_uptake_rate)

if __name__ == "__main__":
    pam_parametrizer = set_up_pamparametrizer(MIN_SUBSTRATE_UPTAKE_RATE, MAX_SUBSTRATE_UPTAKE_RATE)

    pam_parametrizer.run(remove_subruns=True, binned = 'false')
# for running:
# python -m Scripts.Testing.pam_parametrizer_toy_model
