import unittest
import pytest
import os
import pandas as pd
import numpy as np
from scipy.stats import linregress
from pathlib import Path
from typing import Union
from deap.base import Fitness

from Modules.genetic_algorithm_parametrization import GAPOUniform
from Modules.utils.sector_config_functions import change_translational_sector_with_config_dict, get_model_simulations_vs_sector, perform_linear_regression

from Scripts.pam_generation import setup_toy_pam


class GeneticAlgorithmMock(GAPOUniform):
    RESULT_DF_FILE = os.path.join('Scripts', 'i2_parametrization', 'Data', 'toy_model_simulations_ga.csv')

    def __init__(self,
                 substrate_uptake_rates = [0.001,0.091]):
        valid_data_df = pd.read_csv(self.RESULT_DF_FILE).round({'R1_ub':3}) #need to round for correct matching to simulations
        pamodel = setup_toy_pam()


        super().__init__(model=pamodel, # Metabolic model,
        enzymes_to_eval = {'E3':{'reaction':'R3','kcats': {'f': 1}, 'sensitivity':0.5}, #should become 5
                           'E4':{'reaction':'R4','kcats': {'f': 0.5}, 'sensitivity':0.2},#should become 0.1
                           'E5':{'reaction':'R5','kcats': {'f': 0.45}, 'sensitivity':0.1}},#should become 0.25
        fitness_class="Fitfun_params_uniform", # filename (or module) of the fitness function class
        mutation_probability = 0.5, # probability with which an individual (solution) is mutated in a generation
        mutation_rate = 0.5, # probability with which an attribute (e.g. gene) of an individual is mutated
        population_size = 5, # number of individuals (solution) per population
        crossover_probability = 0.8, # probability with which two indivduals/offsprings are crossed over
        number_generations = 3, # number of consecutive generations per gene flow event
        number_gene_flow_events = 2, # number of gene flow events, i.e. merging of multiple
                                     # populations independently evolved on parallel workers
        init_attribute_probability=0.001, # probability with which attributes of initial individuals are mutated
        processes = 1, # number of parallel workers
        time_limit = 600, # time limit in seconds
        filename_save=f"test_ga_parametrization", # filename for saving results after every gene flow event
        folderpath_save=Path(r"Results"), # path for saving results
        overwrite_intermediate_results=True, # if true, saved intermediate results are overwritten
        valid_data ={'R1':valid_data_df},
        sigma_denominator= 10,
        objective_id = 'R7',
        substrate_uptake_rates = {'R1':substrate_uptake_rates},
        substrate_uptake_id = 'R1')

    def get_initial_population(self):
        toolbox = self._init_deap_toolbox()
        population = self.ga.init_pop(toolbox, self.population_size, True)
        return population



def test_genetic_algorithm_initiates_population():
    # Arrange
    sut = GeneticAlgorithmMock()
    # Act
    population = sut.get_initial_population()

    # Assert
    assert sut.population_size == len(population)

def test_genetic_algorithm_finds_best_individual():
    # Arrange
    sut = GeneticAlgorithmMock()
    population = sut.get_initial_population()

    #define the best individual
    elite_actual = population[3]
    elite_actual.fitness.values = [1]
    population[3] = elite_actual

    # Act
    elite_according_to_ga = sut.ga._get_best_individual_from_population(population)

    # Assert
    assert elite_actual == elite_according_to_ga

def test_genetic_algorithm_clones_elite_properly():
    # Arrange
    sut = GeneticAlgorithmMock()
    toolbox = sut._init_deap_toolbox()
    population = sut.get_initial_population()

    # define the best individual
    elite = population[3]

    # Act
    elite_cloned = sut.ga._clone_elite([elite], toolbox)

    # Assert
    assert elite is not elite_cloned[0] # it should be a clone, thus not the same instance!
    assert elite.kcat_list == elite_cloned[0].kcat_list
    assert elite.fitness._wsum() == pytest.approx(elite_cloned[0].fitness._wsum(), abs=1e-3)


def test_genetic_algorithm_adjust_kcat_correctly():
    # Arrange
    sut = GeneticAlgorithmMock()
    population = sut.get_initial_population()
    individual_ut = population[0]
    kcat_old = individual_ut.kcat_list
    kcat_new = [5,0.1,0.25]
    individual_ut.kcat_list = kcat_new
    constraint_names = [f'EC_{enzyme}_f' for enzyme in individual_ut.enzymes_to_eval]
    reaction_names = [f"CE_{rxn_id}_{enz_id}" for rxn_id, enz_id in zip(individual_ut.reactions, individual_ut.enzymes_to_eval)]
    model = sut.FitEval.model

    # Act
    sut.FitEval._change_kcat_values_for_individual(individual_ut)
    #get kcats from model
    kcats_after_ga_adjustment = []
    for rxn_id, constraint_id in zip(reaction_names, constraint_names):
        rxn = model.reactions.get_by_id(rxn_id)
        coeff = model.constraints[constraint_id].get_linear_coefficients([rxn.forward_variable])[rxn.forward_variable]
        kcats_after_ga_adjustment.append(coeff)
        # kcats_after_ga_adjustment += [1/coeff/(3600*1e-6)]#unit conversion

    # Assert
    assert kcat_old != pytest.approx(kcats_after_ga_adjustment,abs=1e-2)
    assert kcat_new == pytest.approx(kcats_after_ga_adjustment,abs=1e-2)

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
    assert new_kcats == [1/ kcat for kcat in individual_to_evaluate.kcat_list]
    assert fitness_validation == pytest.approx(fitness_simulated, abs=1e-6)

def test_genetic_algorithm_toolbox_evaluate_function_gives_right_output():
    # Arrange
    sut = GeneticAlgorithmMock()
    toolbox = sut._init_deap_toolbox()
    population = sut.get_initial_population()
    kcat_lists = [[5,0.1,0.25],[1,0.5,0.45],[0.1,0.1,0.1],[1,1,1], [3,3,3]]
    for individual, kcat_list in zip(population, kcat_lists):
        # individual.kcat_list = [kcat/(3600*1e-6) for kcat in kcat_list]
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
    sut.FitEval.translational_sector_config = {'R1': {'slope':0.01, 'intercept':0.01*1e-3},
                                             'R9': {'slope':0, 'intercept':0.01*1e-3}}

    # Set up a population for comparison
    toolbox = sut._init_deap_toolbox()
    population = toolbox.population(n=3)
    population[0].kcat_list = [1 / kcat for kcat in new_kcats]

    # Act
    fitnesses_from_ga = [fit._wsum() for fit in map(toolbox.evaluate, population)]
    # Expect a perfect fit, as the same model was used to generate the validation data
    assert fitnesses_from_ga[0] == 1

def test_fitness_evaluation_configures_translational_sector_correctly():
    # Arrange
    sut = GeneticAlgorithmMock()
    slope = 1*1e-3
    intercept = 0.01*1e-3
    tps_0 = sut.FitEval.model.sectors.get_by_id('TranslationalProteinSector').intercept
    tot_prot = sut.FitEval.model.constraints[sut.FitEval.model.TOTAL_PROTEIN_CONSTRAINT_ID].ub + tps_0

    # Apply
    change_translational_sector_with_config_dict(sut.FitEval.model,
                                                 {'slope': slope, 'intercept': intercept},
                                                 'R1')
    # Assert
    rxn = sut.FitEval.model.reactions.R1
    coeff = sut.FitEval.model.constraints[sut.FitEval.model.TOTAL_PROTEIN_CONSTRAINT_ID].get_linear_coefficients([rxn.forward_variable])[rxn.forward_variable]
    assert sut.FitEval.model.constraints[sut.FitEval.model.TOTAL_PROTEIN_CONSTRAINT_ID].ub == tot_prot-intercept*1e3
    assert -slope == coeff



##########################################################################################################################
# HELPER FUNCTIONS
##########################################################################################################################

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
                               reference_data_file:Union[str, pd.DataFrame] = 'Scripts/i2_parametrization/Data/toy_model_simulations_ga.csv',
                               substrate_rxn:str = 'R1_ub') -> float:
    """
    Evaluate the fitness of the toymodel compared to the reference dataset generated using kcat_fwd = [1, 0.5, 5, 0.1, 0.25, 1.5]
    :return: float: error average difference of validation and result for the total of substrate uptake range and available reactiosn
    """
    if isinstance(reference_data_file, str):
        validation_data = pd.read_csv(reference_data_file).round({'R1_ub':3})
    else:
        validation_data = reference_data_file
    simulation_results = run_simulations(toy_model, substrate_rates)
    error = []
    for reaction_id in validation_data.columns[2:]:
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
                                            reaction_id='R9', substrate_uptake_rates=[-1e-3, -1e-2]):
    reactions_to_save = [reaction_id, 'R7', 'R8', 'R1']
    fluxes_df = pd.DataFrame(columns = [f'{reaction_id}_ub']+reactions_to_save)

    toy_model.change_reaction_bounds('R1', -1e6, 0)
    for rate in substrate_uptake_rates:
        toy_model.change_reaction_bounds(reaction_id, rate, 0)
        toy_model.optimize()
        flux_results = [rate] + [toy_model.reactions.get_by_id(rxn_id).flux for rxn_id in reactions_to_save]
        fluxes_df.loc[len(fluxes_df)] = flux_results
    toy_model.change_reaction_bounds('R9', 0, 1e6)
    toy_model.change_reaction_bounds('R1', -1e6, 1e6)
    return fluxes_df, reactions_to_save
