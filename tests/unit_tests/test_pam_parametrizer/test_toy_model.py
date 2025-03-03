import pandas as pd

import pytest
from Scripts.i2_parametrization.pam_parametrizer_toy_model import set_up_pamparametrizer
from Scripts.pam_generation import setup_toy_pam
from tests.unit_tests.test_genetic_algorithm_parametrization.test_ga_params import (evaluate_toy_model_fitness,
                                                                                    get_toy_model_simulations_other_csource)
from tests.pam_parametrizer_mock import PAMParametrizerMock
from tests.unit_tests.test_pam_parametrizer.test_pam_parametrizer import save_simulated_fluxes_in_pamparametrizer_for_different_carbon_sources
from tests.unit_tests.test_genetic_algorithm_parametrization.test_ga_params import GeneticAlgorithmMock

FINAL_ENZYMES2KCAT = {'E1':{'CE_R1_E1':{'f': 1/(3600*1e-6), 'b':1/(3600*1e-6)}},
                          'E2': {'CE_R2_E2': {'f': 0.5 / (3600 * 1e-6), 'b': 0.5 / (3600 * 1e-6)}},
                          'E3':{'CE_R3_E3':{'f': 5/(3600*1e-6), 'b':5/(3600*1e-6)}},
                      'E4':{'CE_R4_E4':{'f': 0.1/(3600*1e-6), 'b':0.1/(3600*1e-6)}},
                      'E5':{'CE_R5_E5':{'f': 0.25/(3600*1e-6), 'b':0.25/(3600*1e-6)}},
                          'E6':{'CE_R6_E6':{'f': 1.5/(3600*1e-6), 'b':1.5/(3600*1e-6)}}}

def test_if_toy_model_parameters_in_pam_parametrizer_are_set_correctly():
    # Arrange
    sut = set_up_pamparametrizer(0.001, 0.1, kcat_fwd = [1, 0.5, 5, 0.1, 0.25, 1.5])
    # change_kcat_to_expected_outcome(sut)

    # Assert
    for enzyme, rxn2kcat_dict in FINAL_ENZYMES2KCAT.items():
        kcat_to_validate = sut._pamodel.enzymes.get_by_id(enzyme).rxn2kcat
        assert all(v == kcat_to_validate[k] for k,v in rxn2kcat_dict.items()) and len(rxn2kcat_dict)==len(kcat_to_validate)

def test_if_running_toy_model_in_pam_parametrizer_gives_correct_results():
    # Arrange
    sut = set_up_pamparametrizer(0.001, 0.1, kcat_fwd = [1, 0.5, 5, 0.1, 0.25, 1.5])
    sut.validation_data.get_by_id('R1')._reactions_to_validate = ['R1', 'R7', 'R8', 'R9']
    sut._init_results_objects()
    bin_id = 1
    bin_information = [0.001,0.1, 1e-2]
    # change_kcat_to_expected_outcome(sut)
    reference_flux_data = sut.validation_data.get_by_id('R1').valid_data.iloc[:,1:]

    # Act
    sut.run_pamodel_simulations_in_bin(bin_id = bin_id, bin_information=bin_information, substrate_uptake_reaction='R1')
    flux_data = sut.parametrization_results.flux_results.get_by_id('R1').fluxes_df.drop('bin', axis =1)
    flux_data = flux_data.rename(columns = {'substrate': 'R1_ub'})

    # Assert
    assert pd.testing.assert_frame_equal(reference_flux_data, flux_data,
                                         check_dtype=False,check_exact=False, atol=1e-3) is None


def test_if_simulation_error_for_multiple_carbon_sources_of_parametrizer_is_same_as_for_genetic_algorithm():
    # Arrange
    #setup PAMparametrizer
    sut_param = PAMParametrizerMock()
    # Get the simulations with the byproduct as carbon source
    other_substrate_reaction = 'R9'
    substrate_uptake_rates = [-1e-2, -1e-3]
    toy_pam = setup_toy_pam()
    sut_param._pamodel = toy_pam

    expected_flux_results, reactions_to_validate = get_toy_model_simulations_other_csource(toy_pam,
                                                                                           other_substrate_reaction,
                                                                                           substrate_uptake_rates)
    # Change the validation data object in the parametrization object
    sut_param.add_new_substrate_source(new_substrate_uptake_id=other_substrate_reaction,
                                 validation_data=expected_flux_results,
                                 substrate_range=substrate_uptake_rates,
                                 reactions_to_validate=reactions_to_validate)

    sut_param._init_validation_df(substrate_uptake_ids = ['R1', 'R9'])
    sut_param.validation_data.get_by_id('R1').sampled_valid_data = sut_param.validation_data.get_by_id('R1').valid_data
    sut_param = save_simulated_fluxes_in_pamparametrizer_for_different_carbon_sources(sut_param)

    # set up genetic algorithm
    sut_ga = set_up_mockga_multiple_csources(expected_flux_results, 'R9', reactions_to_validate, substrate_uptake_rates)
    toolbox = sut_ga._init_deap_toolbox()
    toy_ga = sut_ga.ga

    sut_ga.FitEval.substrate_uptake_rates['R1'] = sut_param.validation_data.get_by_id('R1').valid_data['R1_ub'].to_list()
    population = toolbox.population(n=3)
    population[0].kcat_list = [1, 0.5 ,0.45]

    # Act
    #pam parametrizer
    r_squared_param = sut_param.calculate_final_error()
    #genetic algorithm
    population = toy_ga.evaluate_pop(population, toolbox)
    individual_to_evaluate = population[0]
    r_squared_ga = individual_to_evaluate.fitness._wsum()

    # Assert
    assert r_squared_ga == pytest.approx(r_squared_param, abs = 1e-3)


##########################################################################################################################
# HELPER FUNCTIONS
##########################################################################################################################
def change_kcat_to_expected_outcome(pam_parametrizer):
    for enz_id, kcat_dict in FINAL_ENZYMES2KCAT.items():
        pam_parametrizer._pamodel.change_kcat_value(enz_id, kcat_dict)

def set_up_mockga_multiple_csources(expected_flux_results, new_substrate_id,
                                    reactions_to_validate, substrate_uptake_rates):
    sut_ga = GeneticAlgorithmMock()
    fiteval = sut_ga.FitEval
    fiteval.valid_data[new_substrate_id] = expected_flux_results
    fiteval.reactions_with_data[new_substrate_id] = reactions_to_validate
    fiteval.substrate_uptake_rates[new_substrate_id] = substrate_uptake_rates
    return sut_ga
