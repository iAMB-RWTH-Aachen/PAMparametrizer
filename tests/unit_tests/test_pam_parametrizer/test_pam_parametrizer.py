import PAModelpy.PAModel
import pandas as pd
import numpy as np
import os

import pytest
from Modules.PAM_parametrizer import ValidationData, HyperParameters, ParametrizationResults
from Modules.PAM_parametrizer import PAMParametrizer
from Scripts.pam_generation import setup_toy_pam


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



        self.parametrization_results.initiate_result_dfs(reactions_to_validate= ['R1', 'R7', 'R8', 'R9'],
                                                         biomass_reaction= ['R7'])


    def set_up_validation_data_mock(self):
        script_directory = os.path.dirname(os.path.abspath(__file__))
        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(script_directory)))
        DATA_DIR = os.path.join(base_dir, 'Scripts', 'Testing', 'Data')
        RESULT_DF_FILE = os.path.join(DATA_DIR, 'toy_model_simulations_ga.csv')
        valid_data_df = pd.read_csv(RESULT_DF_FILE)

        validation_data = ValidationData(valid_data_df)
        validation_data._reactions_to_plot = ['R1', 'R7', 'R8', 'R9']
        return validation_data

    def set_up_hyperparameter_mock(self):
        hyperparams = HyperParameters
        hyperparams.number_of_kcats_to_mutate = 3
        return hyperparams



##########################################################################################################################
#TEST FUNCTIONS
##########################################################################################################################

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
    start, stop, step = bin_to_split_info[0], bin_to_split_info[1], bin_to_split_info[2]

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

def test_pam_parametrizer_saves_results_of_simulation_correctly():
    #arrange
    sut = PAMParametrizerMock()
    bin_information = {1:[0.001, 0.002, 0.001/5]}
    pamodel = setup_toy_pam()

    #act
    sut.run_pamodel_simulations_in_bin(bin_information)
    fluxes, esc_coeff_df = run_pamodel_binned(pamodel, bin_information)
    #need to make sure the order is correct
    esc_coeff_df = esc_coeff_df[1].reindex(columns=sut.parametrization_results.esc_df.columns)

    #assert
    assert (esc_coeff_df == sut.parametrization_results.esc_df).all().all()

def test_pam_parametrizer_splits_bins_when_large_variability_in_esc():
    #arrange
    sut = PAMParametrizerMock()
    esc_topn_df = pd.DataFrame({'substrate': list(range(0,10,1)), 'E1':list(range(0,10,1)),
                                'E2': [0]*10,'E3': [1]*10})
    bin_id = 1
    #act
    sut.determine_bin_to_split(esc_topn_df, bin_id)
    sut_result = sut.parametrization_results.bins_to_change
    split_result = sut_result[sut_result['bin'] == bin_id]['split'].iloc[0]
    #assert
    assert split_result


def test_pam_parametrizer_does_not_splits_bins_when_little_variability_in_esc():
    # arrange
    sut = PAMParametrizerMock()
    esc_topn_df = pd.DataFrame({'substrate': list(range(0, 10, 1)), 'E1':[5]*10,
                                'E2': [0] * 10, 'E3': [1] * 10})
    bin_id = 1
    # act
    for enzyme_id in esc_topn_df.columns:
        if enzyme_id == 'substrate': continue

    sut.determine_bin_to_split(esc_topn_df, bin_id)
    sut_result = sut.parametrization_results.bins_to_change
    # assert
    assert len(sut_result) == 0


def test_pam_parametrizer_calculates_r_squared_correctly():
    #arange
    sut = PAMParametrizerMock()
    bin_id = 1
    bin_information = {bin_id:[0.001, 0.002, 0.001/5]}
    pamodel = setup_toy_pam()

    #get exactly the simulation results which we would obtain in the model
    fluxes, esc_coeff_df = run_pamodel_binned(pamodel, bin_information)
    validation_data_mock = fluxes[bin_id]
    sut.run_pamodel_simulations_in_bin(bin_information)

    #act
    r_squared_to_validate = sut._calculate_r_squared_for_reaction('R7', validation_data_mock, bin_information[bin_id], bin_id)

    #assert
    #manual calculation results in an R^2 of 0.88 for this bin
    assert 0.88 == pytest.approx(r_squared_to_validate, 1e-1)

def test_pam_parametrizer_determines_most_sensitive_enzymes_correctly():
    # Arrange
    columns_to_be_present = ['enzyme_id', 'mean', 'bin', 'substrate', 'rxn_id']
    enzymes_to_be_selected = ['E2', 'E5', 'E1']

    sut = PAMParametrizerMock()
    bin_id = 1
    bin_information = {bin_id: [0.08, 0.09, 0.01 / 5]}
    #we need to run simulations before we can parse the results
    sut.run_pamodel_simulations_in_bin(bin_information)

    # Act
    esc_topn_df = sut.determine_most_sensitive_enzymes(bin_id,
                                                       nmbr_kcats_to_pick = 3)
    esc_topn_minimal = esc_topn_df.drop_duplicates(subset = 'enzyme_id')

    # Assert
    #check if the columns required for further processing are available in the dataframe
    for column in columns_to_be_present:
        assert column in esc_topn_df.columns
    #check if the correct enzymes are selected
    assert 3 == len(esc_topn_minimal)
    for enzyme in enzymes_to_be_selected:
        assert enzyme in esc_topn_minimal.enzyme_id.to_list()
    #check if the enzymes are in the right order
    assert esc_topn_minimal.absolute_esc_mean.iloc[0] > esc_topn_minimal.absolute_esc_mean.iloc[-1]

def test_pam_parametrizer_parses_enzymes_to_evaluate_correctly():
    # Arrange
    sut = PAMParametrizerMock()
    esc_topn_df_dummy = pd.DataFrame({
        'bin': [1,1,1],
        'enzyme_id': ['E3', 'E4', 'E5'],
        'rxn_id': ['R3', 'R4', 'R5'],
        'mean': [0.5, 0.2, 0.1]
    })

    # unit correction to be consistent with the way it is saved (see toymodel setup function)
    enzymes_to_evaluate_validation = {'E3':{'reaction':'R3','kcat':1/(3600 * 1e-6), 'sensitivity':0.5},
                                      'E4':{'reaction':'R4','kcat':0.5/(3600 * 1e-6), 'sensitivity':0.2},
                                      'E5':{'reaction':'R5','kcat':0.45/(3600 * 1e-6), 'sensitivity':0.1}}

    # Act
    enzymes_to_evaluate_to_test = sut._parse_enzymes_to_evaluate(esc_topn_df_dummy)

    # Assert
    assert all(enzymes_to_evaluate_to_test[key] == value for key, value in enzymes_to_evaluate_validation.items())

def test_if_genetic_algorithm_runs():
    # Arrange
    sut = PAMParametrizerMock()
    esc_topn_df_dummy = pd.DataFrame({
        'bin': [1, 1, 1],
        'enzyme_id': ['E3', 'E4', 'E5'],
        'rxn_id': ['R3', 'R4', 'R5'],
        'mean': [0.5, 0.2, 0.1]
    })
    bin_info = [0.001, 0.002, 0.001/5]
    filename_extension = 'test'
    full_file_path = os.path.join(os.getcwd(),'Results', sut.hyperparameters.genetic_algorithm_filename_base + filename_extension)

    # Act
    sut.run_genetic_algorithm(bin_info, esc_topn_df_dummy, filename_extension)

    # Assert
    # check if filenames are present
    assert all(os.path.exists(full_file_path + file_type) for file_type in ['.json', '.xlsx', '.pickle']),\
        'Genetic algorithm did not create output'
    #remove files
    [os.remove(full_file_path + file_type) for file_type in ['.json', '.xlsx', '.pickle']]

def test_if_pam_parametrizer_gets_all_genetic_algorithm_json_files():
    # Arrange
    sut = PAMParametrizerMock()
    esc_topn_df_dummy = pd.DataFrame({
        'bin': [1, 1, 1],
        'enzyme_id': ['E3', 'E4', 'E5'],
        'rxn_id': ['R3', 'R4', 'R5'],
        'mean': [0.5, 0.2, 0.1]
    })
    bin_information = {1:[0.001, 0.002, 0.001 / 5], 2: [0.001, 0.002, 0.001 / 5]}
    filename_extension = 'test'
    full_file_path = os.path.join(os.getcwd(), 'Results',
                                  sut.hyperparameters.genetic_algorithm_filename_base + filename_extension)
    for bin_id, bin_info in bin_information.items():
        filename_extension_bin = filename_extension + str(bin_id)
        sut.run_genetic_algorithm(bin_info, esc_topn_df_dummy, filename_extension_bin)

    # Act
    json_files = sut._get_genetic_algorithm_json_files()

    # Assert
    assert len(bin_information) == len(json_files)
    # remoce the produced files
    for bin_id in bin_information.keys():
        [os.remove(full_file_path + str(bin_id) +file_type) for file_type in ['.json', '.xlsx', '.pickle']]



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

    return valid_binned_substrate

def run_pamodel_binned(pamodel:PAModelpy.PAModel.PAModel, bin_information:dict):
    fluxes = {}
    esc = {}
    for bin_id, bin_info in bin_information.items():
        esc_df = pd.DataFrame(columns=['bin', 'substrate', 'rxn_id', 'enzyme_id'])
        flux_df = pd.DataFrame(columns= ['substrate','R1', 'R7', 'R8', 'R9'])
        start, stop, step = bin_info[0], bin_info[1], bin_info[2]
        for substrate_uptake_rate in np.arange(start, stop, step):
            pamodel.change_reaction_bounds(rxn_id='R1',
                                                lower_bound=0, upper_bound=substrate_uptake_rate)
            # solve the model
            pamodel.optimize()

            if pamodel.solver.status == 'optimal' and pamodel.objective.value != 0:

                enzyme_sensitivity_coeff = pamodel.enzyme_sensitivity_coefficients
                enzyme_sensitivity_coeff['bin'] = [bin_id] * len(enzyme_sensitivity_coeff)
                enzyme_sensitivity_coeff['substrate'] = [substrate_uptake_rate] * len(enzyme_sensitivity_coeff)
                esc_df = pd.concat([esc_df, enzyme_sensitivity_coeff])

                #save fluxes
                results_row = [substrate_uptake_rate]
                for rxn in ['R1', 'R7', 'R8', 'R9']:
                    results_row += [pamodel.reactions.get_by_id(rxn).flux]

                flux_df.loc[len(flux_df)] = results_row

        esc[bin_id] = esc_df
        fluxes[bin_id] = flux_df

    return fluxes, esc

