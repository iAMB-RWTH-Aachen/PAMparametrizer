import time

import PAModelpy.PAModel
import pandas as pd
import numpy as np
import os
from typing import Callable

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
        DATA_DIR = os.path.join(os.getcwd(), 'Scripts', 'Testing', 'Data')
        RESULT_DF_FILE = os.path.join(DATA_DIR, 'toy_model_simulations_ga.csv')
        valid_data_df = pd.read_csv(RESULT_DF_FILE)

        validation_data = ValidationData(valid_data_df)
        validation_data._reactions_to_plot = ['R1', 'R7', 'R8', 'R9']
        validation_data._reactions_to_validate = ['R1', 'R7', 'R8', 'R9']
        return validation_data

    def set_up_hyperparameter_mock(self):
        hyperparams = HyperParameters
        hyperparams.threshold_iteration = 3
        hyperparams.number_of_kcats_to_mutate = 3
        hyperparams.genetic_algorithm_hyperparams['number_generations'] = 2
        hyperparams.genetic_algorithm_filename_base = 'genetic_algorithm_run_test_'
        hyperparams.genetic_algorithm_hyperparams['print_progress'] = False
        return hyperparams



##########################################################################################################################
#TEST FUNCTIONS
##########################################################################################################################

def test_pam_parametrizer_binning_substrate_uptake_rates_without_splitting():
    # arrange
    sut = PAMParametrizerMock()
    valid_binned_substrate = bin_substrate()

    # act
    sut_binned_substrate = sut._bin_substrate_uptake_rates()
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
    sut_binned_substrate = sut._bin_substrate_uptake_rates()

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
    esc_topn = pd.DataFrame({'substrate': list(range(0,10,1)), 'E1':list(range(0,10,1)),
                                'E2': [0]*10,'E3': [1]*10})

    esc_topn_df = pd.melt(esc_topn, id_vars=['substrate'], var_name='enzyme_id', value_name='coefficient')

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
    esc_topn = pd.DataFrame({'substrate': list(range(0, 10, 1)), 'E1': [5] * 10,
                             'E2': [0] * 10, 'E3': [1] * 10})

    esc_topn_df = pd.melt(esc_topn, id_vars=['substrate'], var_name='enzyme_id', value_name='coefficient')

    bin_id = 1
    # act
    sut.determine_bin_to_split(esc_topn_df, bin_id)
    sut_result = sut.parametrization_results.bins_to_change

    # assert
    assert 0 == len(sut_result)


def test_pam_parametrizer_calculates_r_squared_correctly():
    #arange
    sut = PAMParametrizerMock()
    bin_id = 1
    bin_information = {bin_id:[0.001, 0.002, 0.001/5]}
    start = 0.001
    stop = 0.002
    pamodel = setup_toy_pam()

    #get exactly the simulation results which we would obtain in the model
    fluxes, esc_coeff_df = run_pamodel_binned(pamodel, bin_information)
    validation_data_mock = fluxes[bin_id]
    validation_data_mock['R1_ub'] = validation_data_mock['R1']
    sut.run_pamodel_simulations_in_bin(bin_information)

    flux_df = sut.parametrization_results.fluxes_df
    # check if we want to calculate the error for a single bin
    flux_df = flux_df[flux_df['bin'] == bin_id]

    #act
    r_squared_to_validate = sut._calculate_r_squared_for_reaction(reaction_id ='R7',
                                                                  validation_data = validation_data_mock,
                                                                  substrate_start = start,
                                                                  substrate_stop = stop,
                                                                  fluxes = flux_df)

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
    full_file_path = os.path.join(os.getcwd(),'Results',
                                  sut.hyperparameters.genetic_algorithm_filename_base + filename_extension)

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
    bin_information =  {1: [0.001, 0.002, 0.001 / 5], 2: [0.001, 0.002, 0.001 / 5]}
    sut = run_mock_genetic_algorithm(sut, bin_information)
    full_file_path = os.path.join(os.getcwd(), 'Results',
                                  sut.hyperparameters.genetic_algorithm_filename_base + 'test')

    # Act
    json_files = sut._get_genetic_algorithm_json_files(subset = 'test')

    # Assert
    assert len(bin_information) == len(json_files)
    # remove the produced files
    for bin_id in bin_information.keys():
        [os.remove(full_file_path + str(bin_id) +file_type) for file_type in ['.json', '.xlsx', '.pickle']]


def tests_pam_parametrizer_parses_enzymes_to_evaluate_for_all_bins_correctly():
    # Arrange
    sut = PAMParametrizerMock()
    #run in regions where we know that there is protein limitation and there are thus non-zero esc
    bin_information = {1: [0.07, 0.08, 0.01 / 5], 2: [0.08, 0.09, 0.01 / 5], 3: [0.09, 0.1, 0.01 / 5]}
    sut.run_pamodel_simulations_in_bin(bin_information)
    enzymes_to_evaluate_expected = ['E2', 'E5', 'E1']

    # Act
    enzymes_to_evaluate_test = sut._determine_enzymes_to_evaluate_for_all_bins(nmbr_kcats_to_pick = 3)

    # Assert
    assert all(enzyme_id in enzymes_to_evaluate_test for enzyme_id in enzymes_to_evaluate_expected), "Not all enzyme IDs are present"


def test_if_restart_genetic_algorithm_runs():
    # Arrange
    sut = PAMParametrizerMock()
    #first run simulations to make sure there are results
    bin_information = {1: [0.07, 0.08, 0.01 / 5], 2: [0.08, 0.09, 0.01 / 5], 3: [0.09, 0.1, 0.01 / 5]}
    sut.run_pamodel_simulations_in_bin(bin_information)
    #run genetic algorithm to get output files for restarting
    sut = run_mock_genetic_algorithm(sut, bin_information)

    full_file_path = os.path.join(os.getcwd(), 'Results',
                                  sut.hyperparameters.genetic_algorithm_filename_base + 'final_run_' + str(sut.iteration))
    full_file_path_test = os.path.join(os.getcwd(), 'Results',
                                  sut.hyperparameters.genetic_algorithm_filename_base)
    # Act
    files_to_remove = sut.restart_genetic_algorithm()
    sut._remove_result_files(files_to_remove)
    json_files = sut._get_genetic_algorithm_json_files(subset = 'final_run')

    # Assert
    assert 1 == len(json_files)
    # remove the produced files
    # [os.remove(full_file_path + file_type) for file_type in ['.json', '.xlsx', '.pickle']]
    # for bin_id in bin_information.keys():
    #     [os.remove(full_file_path_test + str(bin_id) +file_type) for file_type in ['.json', '.xlsx', '.pickle']]

def test_pam_parametrizes_reparametrizes_enzymes_correctly():
    # Arrange
    sut = PAMParametrizerMock()
    esc_topn_df_dummy = pd.DataFrame({
        'bin': [1, 1, 1],
        'enzyme_id': ['E3', 'E4', 'E5'],
        'rxn_id': ['R3', 'R4', 'R5'],
        'mean': [0.5, 0.2, 0.1]
    })
    bin_info = [0.001, 0.002, 0.001 / 5]
    filename_extension = f'final_run_{sut.iteration}'
    full_file_path = os.path.join(os.getcwd(), 'Results',
                                  sut.hyperparameters.genetic_algorithm_filename_base + filename_extension)
    sut.run_genetic_algorithm(bin_info, esc_topn_df_dummy, filename_extension)

    best_individual_kcat_df = sut._get_mutated_kcat_values_from_genetic_algorithm()
    kcats_expected = [[row['id'] ,row['rxn_id'], row['value']] for i, row in best_individual_kcat_df.iterrows()]

    # Act
    sut.reparametrize_pam()

    # Assert
    for kcat_info in kcats_expected:
        enzyme_id, rxn_id, kcat_expected = kcat_info[0], kcat_info[1], kcat_info[2]
        kcat_test = 1/sut.pamodel.enzymes.get_by_id(enzyme_id).get_kcat_values([rxn_id])['f']/(3600*1e-6) #adjust for units adjustment in the PAModel
        assert kcat_expected == sut.pamodel.enzymes.get_by_id(enzyme_id).get_kcat_values([rxn_id])['f']

    # remove the produced files
    [os.remove(full_file_path + file_type) for file_type in ['.json', '.xlsx', '.pickle']]

def test_pam_parametrizer_plots_validation_data():
    # Arrange
    sut = PAMParametrizerMock()

    #remove the produced files
    full_file_path = os.path.join(os.getcwd(), 'Results',
                                  sut.hyperparameters.genetic_algorithm_filename_base)

    # Act
    fig, axs = sut.plot_valid_data()


def test_pam_parametrizer_plots_progress():
    # Arrange
    sut = PAMParametrizerMock()
    fig, axs = sut.plot_valid_data()
    # Act
    fig = sut.plot_simulation(fig, axs)

def test_pam_parametrizer_calculates_final_error_correctly():
    # Arrange
    sut = PAMParametrizerMock()

def test_pam_parametrizer_if_diagnostics_are_saved_to_dataframe():
    # Arrange
    bin_id = 1
    bin_info = [0.001, 0.002, 0.001/5]

    sut = PAMParametrizerMock()
    sut._init_results_objects()
    sut.iteration = 1

    start_time = time.perf_counter()
    sut.process_bin(bin_info, bin_id)
    computational_time = time.perf_counter() - start_time

    sut.final_error = 1
    results_filename = (os.path.join(os.getcwd(), 'Results', sut.hyperparameters.genetic_algorithm_filename_base + 'iteration_'+
                                     str(sut.iteration) + '_bin_1.xlsx'))

    full_file_path = os.path.join(os.getcwd(), 'Results',
                                  sut.hyperparameters.genetic_algorithm_filename_base)
    # Act
    sut.save_diagnostics(computational_time,
                         results_filename)

    # Assert
    assert_run_diagnostics_are_saved(sut.parametrization_results,
                                     number_of_best_individuals= sut.hyperparameters.number_of_kcats_to_mutate)

    #remove produced files
    [os.remove(results_filename[:-5] + file_type) for file_type in ['.json', '.xlsx', '.pickle']]

#TODO
def test_pam_parametrizer_runs_full_workflow():
    # Arrange
    sut = PAMParametrizerMock()
    sut.min_substrate_uptake_rate = 0.07
    sut.max_substrate_uptake_rate = 0.09
    sut.hyperparameters.threshold_error = 1

    # Act
    sut.run()
    #remove result files (incl figure)
    filename_extension = f'final_run_{sut.iteration}'
    full_file_path = os.path.join(os.getcwd(), 'Results',
                                  sut.hyperparameters.genetic_algorithm_filename_base + filename_extension)
    [os.remove(full_file_path + file_type) for file_type in ['.json', '.xlsx', '.pickle']]
    os.remove(sut.result_figure_file)

    # Assert
    # if it runs all is fine, only testing functionality
    assert True


###########################################################################
#HELPER FUNCTIONS
###########################################################################
def bin_substrate() -> dict:
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

def run_pamodel_binned(pamodel:PAModelpy.PAModel.PAModel, bin_information:dict) -> tuple:
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

def run_mock_genetic_algorithm(sut: PAMParametrizerMock,
                               bin_information: dict = {1:
                                                             [0.001, 0.002, 0.001 / 5], 2: [0.001, 0.002, 0.001 / 5]
                                                         }) -> PAMParametrizerMock:
    esc_topn_df_dummy = pd.DataFrame({
        'bin': [1, 1, 1],
        'enzyme_id': ['E3', 'E4', 'E5'],
        'rxn_id': ['R3', 'R4', 'R5'],
        'mean': [0.5, 0.2, 0.1]
    })

    filename_extension = 'test'

    for bin_id, bin_info in bin_information.items():
        filename_extension_bin = filename_extension + str(bin_id)
        sut.run_genetic_algorithm(bin_info, esc_topn_df_dummy, filename_extension_bin)

    return sut

def assert_run_diagnostics_are_saved(parametrization_results_object: Callable,
                                     number_of_best_individuals: int) -> None:
    assert len(parametrization_results_object.best_individuals) == number_of_best_individuals
    assert len(parametrization_results_object.computational_time) == 1
    assert len(parametrization_results_object.final_errors) == 1
