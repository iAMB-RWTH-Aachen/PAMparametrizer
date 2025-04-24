import os
import pandas as pd
from cobra import DictList

from PAModelpy.utils import set_up_pam
from Scripts.pam_generation import setup_toy_pam
from Modules.PAM_parametrizer import PAMParametrizer
from Modules.PAM_parametrizer import ValidationData, HyperParameters, ParametrizationResults, SectorConfig


max_substrate_uptake_rate = 0.1
min_substrate_uptake_rate = 0.001

###########################################################################################################################
#MOCK OBJECTS
###########################################################################################################################

class PAMParametrizerMock(PAMParametrizer):
    def __init__(self):
        kcat_fwd = [i/5 for i in [1, 0.5, 1, 0.5, 0.45, 1.5]]
        toy_pam = setup_toy_pam(kcat_fwd)
        toy_pam.ENZYME_ID_REGEX = r'E([0-9]|[1-9][0-9])'
        validation_data = self.set_up_validation_data_mock()
        hyperparameters = self.set_up_hyperparameter_mock()

        super().__init__(pamodel=toy_pam,
                         validation_data=[validation_data],
                         hyperparameters=hyperparameters,
                         substrate_uptake_id = 'R1',
                         max_substrate_uptake_rate=max_substrate_uptake_rate,
                         min_substrate_uptake_rate = min_substrate_uptake_rate)

        self.result_figure_file = os.path.join('Results', '2_parametrization', 'progress','pam_parametrizer_progress_test.png')


        self.parametrization_results.initiate_result_dfs(reactions_to_validate={'R1':['R1', 'R7', 'R8', 'R9']})


    def set_up_validation_data_mock(self):
        DATA_DIR = os.path.join('tests', 'data')
        RESULT_DF_FILE = os.path.join(DATA_DIR, 'toy_model_simulations_ga.csv')
        valid_data_df = pd.read_csv(RESULT_DF_FILE).round({'R1_ub': 3})

        validation_data = ValidationData(valid_data_df,
                                         'R1',
                                         [min_substrate_uptake_rate, max_substrate_uptake_rate])
        validation_data.sampled_valid_data = valid_data_df
        validation_data._reactions_to_plot = ['R1', 'R7', 'R8', 'R9']
        validation_data._reactions_to_validate = ['R1', 'R7', 'R8', 'R9']



        validation_data.sector_configs = {'TranslationalProteinSector':SectorConfig(
            sectorname = 'TranslationalProteinSector',
            slope = 0.01*1e-3,
            intercept = 0.01*1e-3,
            substrate_range = [-1e-3,-2*1e-3]
        )}
        return validation_data

    def set_up_hyperparameter_mock(self):
        hyperparams = HyperParameters
        hyperparams.threshold_iteration = 3
        hyperparams.number_of_kcats_to_mutate = 3
        hyperparams.genetic_algorithm_hyperparams['number_generations'] = 2
        hyperparams.genetic_algorithm_filename_base = 'genetic_algorithm_run_test_'
        hyperparams.genetic_algorithm_hyperparams['print_progress'] = False
        return hyperparams


class PAMParametrizerMockEcoli(PAMParametrizer):
    '''
    Setting up pam parametrizer for the iML1515 PAM with only the translational protein sector as a predefined sector
    '''
    def __init__(self):
        pam = set_up_pam(os.path.join('tests', 'data', 'proteinAllocationModel_iML1515_EnzymaticData_dummy.xlsx'),
                         sensitivity = False)
        validation_data = self.set_up_validation_data(pam)
        hyperparameters = self.set_up_hyperparameter()

        super().__init__(pamodel=pam,
                         validation_data=[validation_data],
                         hyperparameters=hyperparameters,
                         substrate_uptake_id = 'EX_glc__D_e',
                         max_substrate_uptake_rate=-0.1,
                         min_substrate_uptake_rate = -11)

        self.result_figure_file = os.path.join('Results', '2_parametrization', 'progress','pam_parametrizer_progress_test.png')


    def set_up_validation_data(self, model) -> list[ValidationData]:

        valid_data_df = pd.read_excel(os.path.join('Data', 'Ecoli_phenotypes', 'Ecoli_phenotypes_py_rev.xls'),
                                      sheet_name='Yields')
        valid_data_df = valid_data_df.rename(
            columns={'EX_glc__D_e': 'EX_glc__D_e_ub'})
        valid_data_df = valid_data_df[valid_data_df.Strain == 'MG1655']
        validation_data = ValidationData(valid_data_df, 'EX_glc__D_e',
                                         [-11, -0.1])
        validation_data._reactions_to_plot = []
        validation_data._reactions_to_validate = []

        validation_data.sector_configs = {'TranslationalProteinSector': {
                    'slope': model.sectors.get_by_id('TranslationalProteinSector').tps_mu[0],
                    'intercept': model.sectors.get_by_id('TranslationalProteinSector').tps_0[0]
                }}
        return DictList([validation_data])

    def set_up_hyperparameter(self):
        hyperparams = HyperParameters
        hyperparams.threshold_iteration = 3
        hyperparams.number_of_kcats_to_mutate = 3
        hyperparams.genetic_algorithm_hyperparams['number_generations'] = 2
        hyperparams.genetic_algorithm_filename_base = 'genetic_algorithm_run_test_'
        hyperparams.genetic_algorithm_hyperparams['print_progress'] = False
        return hyperparams