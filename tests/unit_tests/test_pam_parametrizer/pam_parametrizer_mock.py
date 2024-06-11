import os
import pandas as pd

from Scripts.pam_generation import setup_toy_pam
from Modules.PAM_parametrizer import PAMParametrizer
from Modules.PAM_parametrizer import ValidationData, HyperParameters, ParametrizationResults


max_substrate_uptake_rate = 0.01
min_substrate_uptake_rate = 0.001

###########################################################################################################################
#MOCK OBJECTS
###########################################################################################################################

class PAMParametrizerMock(PAMParametrizer):
    def __init__(self):
        kcat_fwd = [i/5 for i in [1, 0.5, 1, 0.5, 0.45, 1.5]]
        toy_pam = setup_toy_pam(kcat_fwd)
        validation_data = self.set_up_validation_data_mock()
        hyperparameters = self.set_up_hyperparameter_mock()

        super().__init__(pamodel=toy_pam,
                         validation_data=validation_data,
                         hyperparameters=hyperparameters,
                         substrate_uptake_id = 'R1',
                         max_substrate_uptake_rate=max_substrate_uptake_rate,
                         min_substrate_uptake_rate = min_substrate_uptake_rate)

        self.result_figure_file = os.path.join('Results', 'pam_parametrizer_progress_test.png')


        self.parametrization_results.initiate_result_dfs(reactions_to_validate={'R1':['R1', 'R7', 'R8', 'R9']},
                                                         biomass_reaction= ['R7'])


    def set_up_validation_data_mock(self):
        DATA_DIR = os.path.join('Scripts', 'Testing', 'Data')
        RESULT_DF_FILE = os.path.join(DATA_DIR, 'toy_model_simulations_ga.csv')
        valid_data_df = pd.read_csv(RESULT_DF_FILE).round({'R1_ub': 3})

        validation_data = ValidationData({'R1':valid_data_df})
        validation_data.sampled_valid_data = {'R1':valid_data_df}
        validation_data._reactions_to_plot = ['R1', 'R7', 'R8', 'R9']
        validation_data._reactions_to_validate = {'R1':['R1', 'R7', 'R8', 'R9']}
        return validation_data

    def set_up_hyperparameter_mock(self):
        hyperparams = HyperParameters
        hyperparams.threshold_iteration = 3
        hyperparams.number_of_kcats_to_mutate = 3
        hyperparams.genetic_algorithm_hyperparams['number_generations'] = 2
        hyperparams.genetic_algorithm_filename_base = 'genetic_algorithm_run_test_'
        hyperparams.genetic_algorithm_hyperparams['print_progress'] = False
        return hyperparams

