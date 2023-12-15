import pandas as pd
import os

from Modules.PAM_parametrizer import ValidationData, HyperParameters, ParametrizationResults
from Modules.PAM_parametrizer import PAMParametrizer
from Scripts.pam_generation import setup_toy_pam


max_substrate_uptake_rate = 0.01
min_substrate_uptake_rate = 0.001
class PAMParametrizerMock(PAMParametrizer):
    def __init__(self):
        toy_pam = setup_toy_pam()
        validation_data = self.set_up_validation_data_mock()
        hyperparameters = self.set_up_hyperparameter_mock()

        super().__init__(pamodel=toy_pam,
                 validation_data=validation_data,
                 hyperparameters=hyperparameters,
                 max_substrate_uptake_rate=max_substrate_uptake_rate,
                         min_substrate_uptake_rate = min_substrate_uptake_rate)


        self.parametrization_results.bins_to_change = pd.DataFrame(columns=['bin', 'split', 'merge'])

    def set_up_validation_data_mock(self):
        DATA_DIR = os.path.join(os.getcwd(), 'Scripts', 'Testing', 'Data')
        RESULT_DF_FILE = os.path.join(DATA_DIR, 'toy_model_simulations_ga.csv')
        valid_data_df = pd.read_csv(RESULT_DF_FILE)

        validation_data = ValidationData
        validation_data.valid_data_df = valid_data_df
        validation_data.reactions_to_plot = ['R1', 'R7', 'R8', 'R9']
        return validation_data

    def set_up_hyperparameter_mock(self):
        hyperparams = HyperParameters
        hyperparams.number_of_kcats_to_mutate = 3
        return hyperparams

def test_pam_parametrizer_binning_substrate_uptake_rates_without_splitting():
    # arrange
    sut = PAMParametrizerMock()
    valid_binned_substrate = bin_substrate()

    # act
    sut_binned_substrate = sut.bin_substrate_uptake_rates()
    # assert
    assert valid_binned_substrate == sut_binned_substrate

def test_pam_parametrizer_binning_substrate_uptake_rates_with_splitting():
    #arrange
    bin_nmbr_to_split = 2
    bin_range = (max_substrate_uptake_rate - min_substrate_uptake_rate) / HyperParameters.number_of_bins
    valid_binned_substrate = bin_substrate()
    bin_to_split_info = valid_binned_substrate[2]
    start, stop, step = bin_to_split_info[0], bin_to_split_info[1],bin_to_split_info[2]

    sut = PAMParametrizerMock()
    sut.parametrization_results.bins_to_change['bin'] = [bin_nmbr_to_split]

    #act
    sut_binned_substrate = sut.bin_substrate_uptake_rates()

    del valid_binned_substrate[bin_nmbr_to_split]
    splitted_bins = {bin_nmbr_to_split: [start, start+bin_range * 0.5, step*0.5],
                     bin_nmbr_to_split + 0.1: [start + bin_range * 0.5, stop, step*0.5]}
    valid_binned_substrate = {**valid_binned_substrate, **splitted_bins}

    #assert
    for key,value in valid_binned_substrate.items():
        assert value == sut_binned_substrate[key]

def test_calculate_error():
    pass

def test_determine_most_sensitive_enzymes():
    pass

def test_reparametrize():
    pass

def test_save_diagnostics():
    pass


###########################################################################
#HELPER FUNCTIONS
###########################################################################
def bin_substrate():
    bin_range = (max_substrate_uptake_rate - min_substrate_uptake_rate) / HyperParameters.number_of_bins
    stepsize = bin_range / HyperParameters.bin_resolution
    substrate_start = 0.001
    valid_binned_substrate = {}
    for i in range(HyperParameters.number_of_bins):
        new_bin = {i: [substrate_start, substrate_start + bin_range, stepsize]}
        valid_binned_substrate = {**valid_binned_substrate, **new_bin}
        # update starting concentration for new bin
        substrate_start += bin_range

    # valid_binned_substrate = {0: [0.001, 0.002, 0.0004], 0.1: [0.002, 0.003, 0.0004], 1: [0.003, 0.005, 0.0004],
    #                           2: [0.005, 0.007, 0.0004], 3: [0.007, 0.009000000000000001, 0.0004],
    #                           4: [0.009000000000000001, 0.011000000000000001, 0.0004]}
    return valid_binned_substrate