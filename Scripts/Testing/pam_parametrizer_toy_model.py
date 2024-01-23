import os
import pandas as pd

print(os.getcwd())

from Modules.PAM_parametrizer import ValidationData, HyperParameters, ParametrizationResults
from Modules.PAM_parametrizer import PAMParametrizer
from Scripts.pam_generation import setup_toy_pam

MAX_SUBSTRATE_UPTAKE_RATE = 0.1
MIN_SUBSTRATE_UPTAKE_RATE = 0.001

def set_up_validation_data():
    DATA_DIR = os.path.join(os.getcwd(), 'Scripts', 'Testing', 'Data')
    RESULT_DF_FILE = os.path.join(DATA_DIR, 'toy_model_simulations_ga.csv')
    valid_data_df = pd.read_csv(RESULT_DF_FILE)

    validation_data = ValidationData(valid_data_df)
    validation_data._reactions_to_plot = ['R1', 'R7', 'R8', 'R9']
    return validation_data

def set_up_hyperparameter():
    hyperparams = HyperParameters
    hyperparams.threshold_iteration = 3
    hyperparams.number_of_kcats_to_mutate = 3
    hyperparams.genetic_algorithm_hyperparams['number_generations'] = 2
    hyperparams.genetic_algorithm_filename_base = 'genetic_algorithm_run_toy_'
    hyperparams.genetic_algorithm_hyperparams['print_progress'] = False
    return hyperparams

def set_up_toy_model():
    kcat_fwd = [i / 10 for i in [1, 0.5, 1, 0.5, 0.45, 1.5]]
    toy_pam = setup_toy_pam(kcat_fwd)
    return toy_pam

def set_up_pamparametrizer(min_substrate_uptake_rate:float, max_substrate_uptake_rate: float):
    toy_pam = set_up_toy_model()
    validation_data = set_up_validation_data()
    hyperparameters = set_up_hyperparameter()

    return PAMParametrizer(pamodel=toy_pam,
                     validation_data=validation_data,
                     hyperparameters=hyperparameters,
                     substrate_uptake_id='R1',
                     max_substrate_uptake_rate=max_substrate_uptake_rate,
                     min_substrate_uptake_rate=min_substrate_uptake_rate)

if __name__ == "__main__":
    pam_parametrizer = set_up_pamparametrizer(MIN_SUBSTRATE_UPTAKE_RATE, MAX_SUBSTRATE_UPTAKE_RATE)
    pam_parametrizer.run()
# for running:
# python -m Scripts.Testing.pam_parametrizer_toy_model
