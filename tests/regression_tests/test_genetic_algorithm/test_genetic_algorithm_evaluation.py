import pandas as pd
import numpy as np
import pytest
import os

from typing import Union
from Modules.PAM_parametrizer import SectorConfig
from Scripts.pam_generation import setup_toy_pam
from tests.unit_tests.test_genetic_algorithm_parametrization.test_ga_params import GeneticAlgorithmMock


def test_genetic_algorithm_calculates_individual_correct_fitness():
    """
    Evaluates the evaluate_pop function, but indirectly the toolbox.evaluate function which is then called
    :return:
    """
    # Arrange
    #setup genetic algorithm and population
    sut = GeneticAlgorithmMock()
    toolbox = sut._init_deap_toolbox()
    toy_ga = sut.ga
    population = toolbox.population(n=3)

    new_kcats = [5,0.1,0.25]
    new_kcats_other = [10, 10, 10]
    population[0].kcat_list = [1/kcat for kcat in new_kcats]
    #also check if other individual is not affected by changing the kcat values
    other_individual = population[2]
    other_individual.kcat_list = [1/kcat for kcat in new_kcats_other]

    # adjust for toy pam altered kcat_values to calculate reference fitness
    kcats = [1, 0.5, 5, 0.1, 0.25, 1.5]
    toy_pam = setup_toy_pam(kcat_fwd=kcats)

    # Act
    population = toy_ga.evaluate_pop(population, toolbox)
    individual_to_evaluate = population[0]
    fitness_simulated = individual_to_evaluate.fitness._wsum()
    fitness_other_indiv_simulated = other_individual.fitness._wsum()

    # Assert
    fitness_validation = evaluate_toy_model_fitness(toy_pam, reference_data_file= sut.RESULT_DF_FILE)

    #1e-6 is solver feasibility tolerance
    assert individual_to_evaluate.fitness is not other_individual.fitness
    assert fitness_other_indiv_simulated != fitness_simulated
    assert new_kcats == [1/kcat for kcat in individual_to_evaluate.kcat_list]
    assert fitness_validation == pytest.approx(fitness_simulated, abs=1e-6)

def test_genetic_algorithm_toolbox_evaluate_function_gives_right_output():
    # Arrange
    sut = GeneticAlgorithmMock()
    toolbox = sut._init_deap_toolbox()
    population = sut.get_initial_population()
    kcat_lists = [[5,0.1,0.25],[1,0.5,0.45],[0.1,0.1,0.1],[1,1,1], [3,3,3]]
    for individual, kcat_list in zip(population, kcat_lists):
        individual.kcat_list = [1/kcat for kcat in kcat_list]

    #calculate actual fitnesses (expected result)
    expected_fitnesses = []
    for kcat_list in kcat_lists:
        # adjust for altered kcat_values
        kcats = [1, 0.5] + kcat_list + [1.5]
        toy_pam = setup_toy_pam(kcat_fwd=kcats)
        expected_fitnesses += [evaluate_toy_model_fitness(toy_pam, reference_data_file=sut.RESULT_DF_FILE)]

    # Act
    fitnesses_from_ga = [fit._wsum() for fit in map(toolbox.evaluate, population)]

    # Assert
    assert len(population) == len(fitnesses_from_ga)
    assert expected_fitnesses == pytest.approx(fitnesses_from_ga, abs = 1e-3)

def test_genetic_algorithm_applies_weighing_scheme_correctly_when_r_squared_is_calculated():
    # Arrange
    sut = GeneticAlgorithmMock()
    weights = {'R8': 5, 'R9':0.1}
    toolbox = sut._init_deap_toolbox()
    toy_ga = sut.ga
    population = toolbox.population(n=1)

    # Act
    population = toy_ga.evaluate_pop(population, toolbox)
    fitness_no_weights = population[0].fitness._wsum()
    sut.FitEval.weights = weights
    population = toy_ga.evaluate_pop(population, toolbox)
    fitness_with_weights = population[0].fitness._wsum()

    # Assert
    assert fitness_no_weights != fitness_with_weights

def test_genetic_algorithms_calculates_error_correct_for_multiple_carbon_sources():
    # Arrange
    sut = GeneticAlgorithmMock()
    # Get the simulations with the byproduct as carbon source
    other_substrate_reaction = 'R9'
    substrate_uptake_rates = [-1e-3, -1e-2]
    new_kcats = [5, 0.1, 0.25]
    toy_pam = setup_toy_pam(kcat_fwd=[1, 0.5] + new_kcats + [1.5])
    expected_flux_results, reactions_to_validate = get_toy_model_simulations_other_csource(toy_pam,
                                                                                           other_substrate_reaction,
                                                                                           substrate_uptake_rates)
    # Change the validation data object in the genetic algorithm
    sut.FitEval.valid_data = {**sut.FitEval.valid_data, **{other_substrate_reaction:expected_flux_results}}
    sut.FitEval.substrate_uptake_rates[other_substrate_reaction] = substrate_uptake_rates
    sut.FitEval.reactions_with_data[other_substrate_reaction] = reactions_to_validate
    sut.FitEval.sector_configs_per_substrate = {
        'R1': {'TranslationalProteinSector':SectorConfig(
            sectorname = 'TranslationalProteinSector',
            slope = 0.01*1e-3,
            intercept = 0.01*1e-3,
            substrate_range = [-1e-3,-2*1e-3]
        )},
        'R9': {'TranslationalProteinSector':SectorConfig(
            sectorname = 'TranslationalProteinSector',
            slope = 0,
            intercept = 0.01*1e-3,
            substrate_range = [-1e-3,-2*1e-3]
        )}}

    # Set up a population for comparison
    toolbox = sut._init_deap_toolbox()
    population = toolbox.population(n=3)
    population[0].kcat_list = [1/kcat for kcat in new_kcats]

    # Act
    fitnesses_from_ga = [fit._wsum() for fit in map(toolbox.evaluate, population)]
    # Expect a perfect fit, as the same model was used to generate the validation data
    assert fitnesses_from_ga[0] == pytest.approx(1, abs=1e-3)

def test_genetic_algorithm_evaluate_reset_old_model_kcats():
    """
    Evaluates the evaluate_pop function, but indirectly the toolbox.evaluate function which is then called
    :return:
    """
    # Arrange
    #setup genetic algorithm
    sut = GeneticAlgorithmMock()
    original_pam = sut.model.copy(copy_with_pickle=True)

    #add different directionalities to see if they are reverted properly
    sut.kcat_list = [1,2,3, 5,5,5]#add some random kcat values
    sut.directions += ['b', 'b', 'b']
    sut.enzymes_to_eval += ['E2', 'E3', 'E4']
    sut.rxns += ['R2', 'R3', 'R4']
    sut._init_deap_individual() # reset deap individual representation

    #set population
    toolbox = sut._init_deap_toolbox()
    toy_ga = sut.ga
    population = toolbox.population(n=3)

    # Act
    toy_ga.evaluate_pop(population, toolbox)

    # Assert
    for rxn_id, enz, direction in zip(sut.rxns, sut.enzymes_to_eval, sut.directions):
        fiteval_rxnvar = [sut.FitEval.model.reactions.get_by_id(rxn_id).forward_variable
                          if direction=='f'
                          else sut.FitEval.model.reactions.get_by_id(rxn_id).reverse_variable
                          ]
        fiteval_coeff = sut.FitEval.model.constraints[
            f'EC_{enz}_{direction}'
        ].get_linear_coefficients(fiteval_rxnvar)[fiteval_rxnvar[0]]

        ori_rxnvar = [original_pam.reactions.get_by_id(rxn_id).forward_variable
                          if direction == 'f'
                          else original_pam.reactions.get_by_id(rxn_id).reverse_variable
                          ]
        ori_coeff = sut.FitEval.model.constraints[
            f'EC_{enz}_{direction}'
        ].get_linear_coefficients(ori_rxnvar)[ori_rxnvar[0]]

        assert ori_coeff == pytest.approx(fiteval_coeff, rel=1e-3)



##################################################################################################################
#HELPER FUNCTIONS
##################################################################################################################

def run_simulations(pamodel, substrate_rates):
    result_df = pd.DataFrame(columns= ['R1_ub','R1', 'R7', 'R8', 'R9'])

    for substrate in substrate_rates:
        pamodel.change_reaction_bounds(rxn_id='R1',
                                       lower_bound=0, upper_bound=substrate)
        print('Running simulations with ', substrate, 'mmol/g_cdw/h of substrate going into the system')
        pamodel.optimize()
        if pamodel.solver.status == 'optimal' and pamodel.objective.value>0:
            results_row = []
            for rxn in ['R1', 'R7', 'R8', 'R9']:
                results_row += [pamodel.reactions.get_by_id(rxn).flux]

            result_df.loc[len(result_df)] = [substrate] + results_row
    return result_df

def evaluate_toy_model_fitness(toy_model, substrate_rates = [0.001, 0.091],
                               reference_data_file:Union[str, pd.DataFrame] = os.path.join(
                                   'tests', 'data', 'toy_model_simulations_ga.csv'
                               ),
                               substrate_rxn:str = 'R1_ub') -> float:
    """
    Evaluate the fitness of the toymodel compared to the reference dataset generated using kcat_fwd = [1, 0.5, 5, 0.1, 0.25, 1.5]
    :return: float: error average difference of validation and result for the total of substrate uptake range and available reactiosn
    """
    if isinstance(reference_data_file, str):
        validation_data = pd.read_csv(reference_data_file).round({'R1_ub':4})
    else:
        validation_data = reference_data_file
    simulation_results = run_simulations(toy_model, substrate_rates)
    error = []
    for reaction_id in validation_data.columns[1:]:
        # Take the absolute value of substrate uptake to avoid issues with reaction directionality
        validation_data[substrate_rxn] = [abs(flux) for flux in validation_data[substrate_rxn]]
        simulated_data = pd.DataFrame({substrate_rxn: [abs(flux) for flux in simulation_results['R1_ub']],
                                       'simulation': simulation_results[reaction_id]})
        ref_data_rxn = pd.merge(validation_data, simulated_data, on=substrate_rxn, how='inner')
        # error: squared difference
        ref_data_rxn = ref_data_rxn.assign(error=lambda x: (x[reaction_id] - x['simulation']) ** 2)

        # calculate R^2:
        data_average = np.nanmean(validation_data[reaction_id])
        residual_ss = np.nansum(ref_data_rxn.error)
        total_ss = np.nansum([(data - data_average) ** 2 for data in ref_data_rxn[reaction_id]])
        # calculating r_squared is only feasible of the numerator and the denomenator are both nonzero
        if (residual_ss == 0) | (total_ss == 0):
            r_squared = 0
        else:
            r_squared = 1 - residual_ss / total_ss

        error += [r_squared]
    return sum(error)/len(error)


def get_toy_model_simulations_other_csource(toy_model,
                                            reaction_id='R9',
                                            substrate_uptake_rates=[-1e-3, -1e-2]):
    reactions_to_save = [reaction_id, 'R7', 'R8', 'R1']
    fluxes_df = pd.DataFrame(columns = [f'{reaction_id}_ub']+reactions_to_save)

    toy_model.change_reaction_bounds('R1', -1e6, 0)
    toy_model.change_sector_parameters(toy_model.sectors.TranslationalProteinSector,
                                       slope=0,
                                       intercept = 0.01*1e-3,
                                       lin_rxn_id=reaction_id)
    for rate in substrate_uptake_rates:
        toy_model.change_reaction_bounds(reaction_id, rate, 0)
        toy_model.optimize()
        flux_results = [rate] + [toy_model.reactions.get_by_id(rxn_id).flux for rxn_id in reactions_to_save]
        fluxes_df.loc[len(fluxes_df)] = flux_results
    toy_model.change_reaction_bounds('R9', 0, 1e6)
    toy_model.change_reaction_bounds('R1', -1e6, 1e6)
    return fluxes_df, reactions_to_save
