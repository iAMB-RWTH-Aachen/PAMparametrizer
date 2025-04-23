import time

import PAModelpy.PAModel
import pandas as pd
import numpy as np
import os
from typing import Callable

import pytest
from Modules.PAM_parametrizer import ValidationData, HyperParameters, ParametrizationResults
from Modules.utils import calculate_r_squared_for_reaction
from Scripts.pam_generation import setup_toy_pam
from tests.regression_tests.test_genetic_algorithm.test_genetic_algorithm_evaluation import (evaluate_toy_model_fitness,
                                                                                    get_toy_model_simulations_other_csource)
from tests.pam_parametrizer_mock import PAMParametrizerMock



max_substrate_uptake_rate = 0.1
min_substrate_uptake_rate = 0.001


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
    assert valid_binned_substrate == sut_binned_substrate['R1']

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
    sut_binned_substrate = sut._bin_substrate_uptake_rates()['R1']

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
    bin_information = [0.001, 0.002, 0.001/5]
    bin_id = 1
    pamodel = setup_toy_pam()

    #act
    sut.run_pamodel_simulations_in_bin('R1', bin_id, bin_information)
    fluxes, esc_coeff_df = run_pamodel_binned(pamodel, {bin_id:bin_information})
    #need to make sure the order is correct
    esc_coeff_df = esc_coeff_df[1].reindex(columns=sut.parametrization_results.esc_df.columns)

    #assert
    assert (esc_coeff_df == sut.parametrization_results.esc_df.round({'bin':0})).all().all()

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
    bin_information = [0.001, 0.002, 0.001/5]
    pamodel = setup_toy_pam()

    #get exactly the simulation results which we would obtain in the model
    fluxes, esc_coeff_df = run_pamodel_binned(pamodel, {bin_id:bin_information})
    validation_data_mock = fluxes[bin_id]
    validation_data_mock['R1_ub'] = validation_data_mock['R1']
    sut.run_pamodel_simulations_in_bin('R1', bin_id, bin_information)

    flux_df = sut.parametrization_results.flux_results.get_by_id('R1').fluxes_df
    # check if we want to calculate the error for a single bin
    flux_df = flux_df[flux_df['bin'] == bin_id]

    #act
    r_squared_to_validate = calculate_r_squared_for_reaction(reaction_id ='R7',
                                                                  validation_data = validation_data_mock.round({'R1_ub': 4}),
                                                                  substrate_uptake_id = 'R1',
                                                                  fluxes = flux_df.round({'R1_ub': 4}))

    #assert
    #manual calculation results in an R^2 of 0.88 for this bin
    assert 1 == pytest.approx(r_squared_to_validate, 1e-1)

def test_pam_parametrizer_determines_most_sensitive_enzymes_correctly():
    # Arrange
    columns_to_be_present = ['enzyme_id', 'mean', 'bin', 'substrate', 'rxn_id']
    enzymes_to_be_selected = ['E2', 'E5', 'E1']

    sut = PAMParametrizerMock()
    bin_id = 1
    bin_information = [0.08, 0.09, 0.01 / 5]
    #we need to run simulations before we can parse the results
    sut.run_pamodel_simulations_in_bin('R1', bin_id, bin_information)

    # Act
    esc_topn_df = sut.determine_most_sensitive_enzymes(bin_id,
                                                       nmbr_kcats_to_pick = 3)
    esc_topn_minimal = esc_topn_df.drop_duplicates(subset = 'enzyme_id')

    # Assert
    #check if the columns required for further processing are available in the dataframe
    for column in columns_to_be_present:
        assert column in esc_topn_df.columns
    #check if the correct enzymes are selected
    print(esc_topn_df)
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
        'rxn_id': ['CE_R3_E3', 'CE_R4_E4', 'CE_R5_E5'],
        'mean': [0.5, 0.2, 0.1]
    })

    enzymes_to_evaluate_validation = {'E3':[{'reaction':'CE_R3_E3','kcats':{'f':1,'b':1}, 'sensitivity':0.5}],
                                      'E4':[{'reaction':'CE_R4_E4','kcats': {'f':0.5, 'b':0.5}, 'sensitivity':0.2}],
                                      'E5':[{'reaction':'CE_R5_E5','kcats':{'f':0.45, 'b':0.45}, 'sensitivity':0.1}]}

    # Act
    enzymes_to_evaluate_to_test = sut._parse_enzymes_to_evaluate(esc_topn_df_dummy)
    # Assert
    for enzid, rxn_list in enzymes_to_evaluate_validation.items():
        rxn_info_validation = rxn_list[0]
        rxn_info_test = enzymes_to_evaluate_to_test[enzid][0]
        #make sure each gpr is only present once
        assert len(rxn_list) == len(enzymes_to_evaluate_to_test[enzid])
        for direction in ['f', 'b']:
            assert rxn_info_test['kcats'][direction] == pytest.approx(1/rxn_info_validation['kcats'][direction], abs=1e-3)


def test_if_genetic_algorithm_runs():
    # Arrange
    sut = PAMParametrizerMock()
    enzymes_to_evaluate= {'E3':[{'reaction':'CE_R3_E3','kcats':{'f':1,'b':1}, 'sensitivity':0.5}],
                                      'E4':[{'reaction':'CE_R4_E4','kcats': {'f':0.5, 'b':0.5}, 'sensitivity':0.2}],
                                      'E5':[{'reaction':'CE_R5_E5','kcats':{'f':0.45, 'b':0.45}, 'sensitivity':0.1}]}

    filename_extension = 'test'
    full_file_path = os.path.join('Results', '2_parametrization',
                                  sut.hyperparameters.genetic_algorithm_filename_base + filename_extension)

    # Act
    sut.run_genetic_algorithm(enzymes_to_evaluate, filename_extension)

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
    full_file_path = os.path.join('Results', '2_parametrization',
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
    for bin_id, bin_info in bin_information.items():
        sut.run_pamodel_simulations_in_bin('R1',bin_id, bin_info)
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
    for bin_id, bin_info in bin_information.items():
        sut.run_pamodel_simulations_in_bin('R1',bin_id, bin_info)
    #run genetic algorithm to get output files for restarting
    sut = run_mock_genetic_algorithm(sut, bin_information)

    # Act
    files_to_remove = sut.restart_genetic_algorithm()
    json_files = sut._get_genetic_algorithm_json_files(subset = 'final_run')
    sut._remove_result_files(files_to_remove)

    # Assert
    assert 1 == len(json_files)

def test_pam_parametrizes_reparametrizes_enzymes_correctly():
    # Arrange
    sut = PAMParametrizerMock()
    enzymes_to_evaluate = {'E3': [{'reaction': 'CE_R3_E3', 'kcats': {'f': 1, 'b': 1}, 'sensitivity': 0.5}],
                                      'E4': [{'reaction': 'CE_R4_E4', 'kcats': {'f': 0.5, 'b': 0.5},
                                             'sensitivity': 0.2}],
                                      'E5': [{'reaction': 'CE_R5_E5', 'kcats': {'f': 0.45, 'b': 0.45},
                                             'sensitivity': 0.1}]}
    filename_extension = f'final_run_{sut.iteration}'
    full_file_path = os.path.join('Results', '2_parametrization',
                                  sut.hyperparameters.genetic_algorithm_filename_base + filename_extension)
    sut.run_genetic_algorithm(enzymes_to_evaluate, filename_extension)

    best_individual_kcat_df, error = sut._get_mutated_kcat_values_from_genetic_algorithm()
    kcats_expected = [[row['id'], row['direction'], row['rxn_id'], row['value']] for i, row in best_individual_kcat_df.iterrows()]
    # Act
    sut.reparametrize_pam()
    # Assert
    for kcat_info in kcats_expected:
        enzyme_id, direction, rxn_id, kcat_expected = kcat_info[0], kcat_info[1], kcat_info[2], kcat_info[3]
        model_kcat_dict = sut._pamodel.enzymes.get_by_id(enzyme_id).get_kcat_values([rxn_id.split("_")[1]])
        if direction in model_kcat_dict.keys():
            kcat_test = (model_kcat_dict[direction]*3600*1e-6) #model units: 1/h *1e6 unit correction, returned ga units: 1/s
            assert 1/kcat_expected == pytest.approx(kcat_test, abs=1e-6)
    # remove the produced files
    [os.remove(full_file_path + file_type) for file_type in ['.json', '.xlsx', '.pickle']]

def test_pam_parametrizer_changes_kcats_same_way_as_genetic_algorithm():
    # Arrange
    kcat = 6
    enzyme_id = 'E1'
    reaction_id = 'R1'
    kcat_dict = {reaction_id: {'f': kcat}}
    enzymes_to_evaluate = {enzyme_id: [{
                        'reaction': reaction_id,
                        'kcats': {'f':1/(kcat* 3600 * 1e-6)},
                        'sensitivity': 0.1}]}
    substrate_rates = {'R1':[1e-3, 1e-2]}

    #set up parametrizer and genetic algorithm objects
    sut = PAMParametrizerMock()
    ga = sut._init_genetic_algorithm(substrate_rates,
                                     enzymes_to_evaluate,
                                     sector_configs_per_substrate={
                                         'R1':sut.validation_data.R1.sector_configs
                                     },
                                     filename_extension='')
    #for the genetic algorithm we need a dummy individual to change the kcat
    toolbox = ga._init_deap_toolbox()
    population = ga.ga.init_pop(toolbox, ga.population_size, True)
    individual = population[0]

    # Act
    sut._change_kcat_value_for_enzyme(enzyme_id='E1',
                                      kcat_dict=kcat_dict)
    kcat_model_sut = get_kcat_values_from_model(sut._pamodel,
                                                enzyme_ids= [enzyme_id],
                                                reaction_names= [f"CE_{reaction_id}_{enzyme_id}"])[0]

    ga.FitEval._change_kcat_values_for_individual(individual)
    kcat_model_ga = get_kcat_values_from_model(ga.FitEval.model,
                                               enzyme_ids= [enzyme_id],
                                               reaction_names= [f"CE_{reaction_id}_{enzyme_id}"])[0]

    # Assert
    assert kcat_model_sut == pytest.approx(kcat_model_ga, abs = 1e-3)
    assert kcat == pytest.approx(kcat_model_sut, abs = 1e-3)
    assert kcat == pytest.approx(kcat_model_ga, abs = 1e-3)



def test_pam_parametrizer_calculates_final_error_correctly():
    # Arrange
    sut = PAMParametrizerMock()
    # Build the PAModel with expected outcome
    sut._pamodel = setup_toy_pam(kcat_fwd = [1, 0.5, 5, 0.1, 0.25, 1.5])
    # get the expected flux distribution
    fluxes, substrate_range = sut.run_simulations_to_plot(substrate_uptake_id='R1')
    # save the flux results
    for simulation_result, substrate_rate in zip(fluxes, substrate_range.keys()):
        sut.parametrization_results.add_fluxes_from_fluxdict(flux_dict=simulation_result,
                                                              bin_id='final',
                                                              substrate_reaction_id='R1',
                                                              substrate_uptake_rate=substrate_rate,
                                                              fluxes_abs=False)

    # Act
    final_error_sut = sut.calculate_final_error()


    #assert
    #R^2 should be 1, because the same model has been used to generate the validation dataset
    assert 1 == pytest.approx(final_error_sut, 1e-2)


def test_pam_parametrizer_if_diagnostics_are_saved_to_dataframe():
    # Arrange
    bin_id = 1
    bin_info = [0.001, 0.002, 0.001/5]

    sut = PAMParametrizerMock()
    sut._init_results_objects()
    sut.iteration = 1

    start_time = time.perf_counter()
    sut.process_bin(bin_id, bin_information=bin_info, substrate_uptake_id='R1')
    computational_time = time.perf_counter() - start_time

    sut.final_error = 1
    results_filename = (os.path.join('Results', '2_parametrization',
                                     sut.hyperparameters.genetic_algorithm_filename_base + '_iteration_'+
                                     str(sut.iteration) + '_bin_1.xlsx'))

    # Act
    sut.save_diagnostics(computational_time,
                         results_filename)

    # Assert
    assert_run_diagnostics_are_saved(sut.parametrization_results,
                                     number_of_best_individuals= sut.hyperparameters.number_of_kcats_to_mutate)

    #remove produced files
    [os.remove(results_filename[:-5] + file_type) for file_type in ['.json', '.xlsx', '.pickle']]


def test_pam_parameterizer_gets_correct_error_for_multiple_carbon_sources():
    # Arrange
    sut = PAMParametrizerMock()
    # Get the simulations with the byproduct as carbon source
    other_substrate_reaction = 'R9'
    substrate_uptake_rates = [-1e-3, -1e-2]
    toy_pam = setup_toy_pam(kcat_fwd=[1, 0.5, 5, 0.1, 0.25, 1.5]) #this was the model used to generate the validation data
    sut._pamodel = toy_pam
    expected_flux_results, reactions_to_validate = get_toy_model_simulations_other_csource(toy_pam,
                                                                                           other_substrate_reaction,
                                                                                           substrate_uptake_rates)

    # Change the validation data object in the parametrization object
    sut.add_new_substrate_source(new_substrate_uptake_id = other_substrate_reaction,
                                 validation_data = expected_flux_results,
                                 substrate_range=substrate_uptake_rates,
                                 reactions_to_validate = reactions_to_validate)
    sut.validation_data.get_by_id(other_substrate_reaction).sampled_valid_data = expected_flux_results

    #get flux data for different carbon sources
    sut = save_simulated_fluxes_in_pamparametrizer_for_different_carbon_sources(sut)

    # Act
    error = sut.calculate_final_error()

    # Assert
    # Expect a perfect fit, as the same model was used to generate the validation data
    assert error == pytest.approx(1, abs = 1e-3)

def test_if_parametrizer_convergence_with_similar_error():
    # Arrange
    sut = PAMParametrizerMock()
    sut.parametrization_results.final_errors = pd.DataFrame({'run_id': [0,1,2,3,4,5],
                                                            'r_squared': [0.11,0.12,0.13,0.13,0.13,0.13]})
    # Act
    converged = sut._error_is_converging()

    # Assert
    assert converged

def test_if_parametrizer_runs_iteration_with_random_enzymes_to_evaluate():
    # Arrange
    sut = PAMParametrizerMock()
    sut._pamodel.objective = "R7"
    sut._pamodel.change_reaction_bounds("R1", 1e-2, 2e-2)
    sut._pamodel.optimize()

    # Act
    files_to_remove = sut.perform_iteration_without_bins(random = True)

    # Assert
    assert len(files_to_remove)>=1
    [os.remove(file) for file in files_to_remove]

#########################################################################################################################
#HELPER FUNCTIONS
#########################################################################################################################
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

def run_pamodel_binned(pamodel:PAModelpy.PAModel, bin_information:dict) -> tuple:
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
                               bin_information: dict = {1:[0.001, 0.002, 0.001 / 5],
                                                        2: [0.001, 0.002, 0.001 / 5]}) -> PAMParametrizerMock:
    enzymes_to_evaluate = {'E3': [{'reaction': 'CE_R3_E3', 'kcats': {'f': 1, 'b': 1}, 'sensitivity': 0.5}],
                                      'E4': [{'reaction': 'CE_R4_E4', 'kcats': {'f': 1 / 0.5, 'b': 1 / 0.5},
                                             'sensitivity': 0.2}],
                                      'E5': [{'reaction': 'CE_R5_E5', 'kcats': {'f': 1 / 0.45, 'b': 1 / 0.45},
                                             'sensitivity': 0.1}]}

    filename_extension = 'test'

    for bin_id, bin_info in bin_information.items():
        filename_extension_bin = filename_extension + str(bin_id)
        sut.run_genetic_algorithm(enzymes_to_evaluate,
                                  filename_extension_bin)

    return sut

def get_kcat_values_from_model(model:PAModelpy.PAModel, enzyme_ids:list, reaction_names: list) -> list:
    # get kcats from model
    constraint_names = [f'EC_{enz_id}_f' for enz_id in enzyme_ids]
    kcats_after_ga_adjustment = []
    for rxn_id, constraint_id in zip(reaction_names, constraint_names):
        rxn = model.reactions.get_by_id(rxn_id)
        coeff = model.constraints[constraint_id].get_linear_coefficients([rxn.forward_variable])[rxn.forward_variable]
        kcats_after_ga_adjustment += [1 / coeff / (3600 * 1e-6)]  # unit conversion
    return kcats_after_ga_adjustment

def assert_run_diagnostics_are_saved(parametrization_results_object: Callable,
                                     number_of_best_individuals: int) -> None:

    assert len(parametrization_results_object.best_individuals) == number_of_best_individuals*2 #times 2 because of f and b for each enzyme
    assert len(parametrization_results_object.computational_time) == 1
    assert len(parametrization_results_object.final_errors) == 1

def save_simulation_results(parametrizer:PAMParametrizerMock,
                            substrate_rates:list, fluxes:list, substrate_uptake_id: str) -> None:
    for simulation_result, substrate_rate in zip(fluxes, substrate_rates):
        parametrizer.parametrization_results.add_fluxes_from_fluxdict(flux_dict=simulation_result,
                                                              bin_id='final',
                                                              substrate_reaction_id=substrate_uptake_id,
                                                              substrate_uptake_rate=substrate_rate,
                                                              fluxes_abs=False)
    return parametrizer

def save_simulated_fluxes_in_pamparametrizer_for_different_carbon_sources(parametrizer):
    for substr_uptake_id in parametrizer.substrate_uptake_ids:
        substr_uptake_rates = parametrizer.validation_data.get_by_id(
            substr_uptake_id).sampled_valid_data[substr_uptake_id+'_ub'].to_list()
        fluxes, substrate_range = parametrizer.run_simulations_to_plot(substrate_uptake_id=substr_uptake_id,
                                                                       substrate_rates=substr_uptake_rates)
        parametrizer = save_simulation_results(parametrizer, substrate_range, fluxes, substr_uptake_id)
    return parametrizer