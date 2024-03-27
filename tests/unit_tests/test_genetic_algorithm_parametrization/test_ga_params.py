import unittest
import pytest
import os
import pandas as pd
from scipy.stats import linregress
from pathlib import Path
from typing import Union
from deap.base import Fitness

from Modules.genetic_algorithm_parametrization import GAPOUniform
from Scripts.pam_generation import setup_toy_pam


class GeneticAlgorithmMock(GAPOUniform):
    RESULT_DF_FILE = os.path.join('Scripts', 'Testing', 'Data', 'toy_model_simulations_ga.csv')

    def __init__(self,
                 substrate_uptake_rates = [0.001,0.091]):
        valid_data_df = pd.read_csv(self.RESULT_DF_FILE)
        pamodel = setup_toy_pam()


        super().__init__(model=pamodel, # Metabolic model,
        enzymes_to_eval = {'E3':{'reaction':'R3','kcat':1, 'sensitivity':0.5}, #should become 5
                           'E4':{'reaction':'R4','kcat':0.5, 'sensitivity':0.2},#should become 0.1
                           'E5':{'reaction':'R5','kcat':0.45, 'sensitivity':0.1}},#should become 0.25
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
        valid_df = valid_data_df,
        sigma_denominator= 10,
        objective_id = 'R7',
        substrate_uptake_rates = substrate_uptake_rates,
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
    reaction_names = individual_ut.reactions
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
    # new_kcats = [kcat/(3600*1e-6) for kcat in [5,0.1,0.25]]
    # new_kcats_other = [kcat / (3600 * 1e-6) for kcat in [10, 10, 10]]

    new_kcats = [5,0.1,0.25]
    new_kcats_other = [10, 10, 10]
    population[0].kcat_list = [1/kcat for kcat in new_kcats]
    #also check if other individual is not affected by changing the kcat values
    other_individual = population[2]
    other_individual.kcat_list = [1/kcat for kcat in new_kcats_other]

    # Act
    population = toy_ga.evaluate_pop(population, toolbox)
    individual_to_evaluate = population[0]
    fitness_simulated = individual_to_evaluate.fitness._wsum()
    fitness_other_indiv_simulated = other_individual.fitness._wsum()

    #adjust for altered kcat_values
    kcats = [1, 0.5,5,0.1,0.25,1.5]
    toy_pam = setup_toy_pam(kcat_fwd = kcats)
    # for key, value in toy_pam.constraints.items():
    #     print(value)
    #
    # for key, value in sut.FitEval.model.constraints.items():
    #     print(value)

    # Assert
    fitness_validation = evaluate_toy_model_fitness(toy_pam, reference_data_file_path = sut.RESULT_DF_FILE)

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
        expected_fitnesses += [evaluate_toy_model_fitness(toy_pam, reference_data_file_path=sut.RESULT_DF_FILE)]

    # Act
    fitnesses_from_ga = [fit._wsum() for fit in map(toolbox.evaluate, population)]

    # Assert
    assert len(population) == len(fitnesses_from_ga)
    assert expected_fitnesses == pytest.approx(fitnesses_from_ga, abs = 1e-3)



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
                               reference_data_file_path:Union[str, pd.DataFrame] = 'Scripts/Testing/Data/toy_model_simulations_ga.csv',
                               substrate_rxn:str = 'R1_ub') -> float:
    """
    Evaluate the fitness of the toymodel compared to the reference dataset generated using kcat_fwd = [1, 0.5, 5, 0.1, 0.25, 1.5]
    :return: float: error average difference of validation and result for the total of substrate uptake range and available reactiosn
    """
    if isinstance(reference_data_file_path, str):
        validation_results = pd.read_csv(reference_data_file_path)
    else:
        validation_results = reference_data_file_path
    simulation_results = run_simulations(toy_model, substrate_rates)
    error = []
    for rxn in validation_results.columns[2:]:
        line = linregress(x=[abs(substrate_rates[0]), abs(substrate_rates[-1])],
                              y=[simulation_results[rxn].iloc[0], simulation_results[rxn].iloc[-1]])
        ref_data_rxn = validation_results.assign(
            simulation=lambda x: line.intercept + line.slope * x[substrate_rxn])
        # simulation mean
        data_average = ref_data_rxn[rxn].mean()
        # error: squared difference
        ref_data_rxn = ref_data_rxn.assign(error=lambda x: (x[rxn] - x['simulation']) ** 2)

        # calculate R^2:
        residual_ss = sum(ref_data_rxn.error)
        total_ss = sum([(data - data_average) ** 2 for data in ref_data_rxn[rxn]])
        # calculating r_squared is only feasible of the numerator and the denomenator are both nonzero
        if (residual_ss == 0) | (total_ss == 0):
            r_squared = 0
        else:
            r_squared = 1 - residual_ss / total_ss

        error += [r_squared]
    return sum(error)/len(error)