import pytest
import os
import pandas as pd
from pathlib import Path

from Modules.PAM_parametrizer import SectorConfig
from Modules.genetic_algorithm_parametrization import GAPOUniform
from Modules.utils.sector_config_functions import (change_proteinsector_relation_from_growth_to_substrate_uptake,
                                                   get_model_simulations_vs_sector,
                                                   perform_linear_regression)

from Scripts.pam_generation import setup_toy_pam


class GeneticAlgorithmMock(GAPOUniform):
    RESULT_DF_FILE = os.path.join('tests', 'data', 'toy_model_simulations_ga.csv')

    def __init__(self,
                 substrate_uptake_rates = [0.001,0.091]
                 ):
        valid_data_df = pd.read_csv(self.RESULT_DF_FILE).round({'R1_ub':3}) #need to round for correct matching to simulations
        pamodel = setup_toy_pam()


        super().__init__(
            model=pamodel, # Metabolic model,
            enzymes_to_eval = {'E3':[{'reaction':'R3','kcats': {'f': 1}, 'sensitivity':0.5}], #should become 5
                               'E4':[{'reaction':'R4','kcats': {'f': 0.5}, 'sensitivity':0.2}],#should become 0.1
                               'E5':[{'reaction':'R5','kcats': {'f': 0.45}, 'sensitivity':0.1}]#should become 0.25
                               },
            sector_configs_per_substrate = {
            'R1': {'TranslationalProteinSector':SectorConfig(
            sectorname = 'TranslationalProteinSector',
            slope = 0.01*1e-3,
            intercept = 0.01*1e-3,
            substrate_range = [-1e-3,-2*1e-3]
            )},
            'R9': {
                'TranslationalProteinSector':SectorConfig(
                    sectorname = 'TranslationalProteinSector',
                    slope = 0,
                    intercept = 0.01*1e-3,
                    substrate_range = [-1e-3,-2*1e-3]
                )}
        },
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

def test_fitness_evaluation_configures_translational_sector_correctly():
    # Arrange
    sut = GeneticAlgorithmMock()
    slope = 1*1e-3
    intercept = 0.01*1e-3
    tps_0 = sut.FitEval.model.sectors.get_by_id('TranslationalProteinSector').intercept
    tot_prot = sut.FitEval.model.constraints[sut.FitEval.model.TOTAL_PROTEIN_CONSTRAINT_ID].ub + tps_0

    # Apply
    change_proteinsector_relation_from_growth_to_substrate_uptake(sut.FitEval.model,
                                                 {'slope': slope, 'intercept': intercept},
                                                 'R1')
    # Assert
    rxn = sut.FitEval.model.reactions.R1
    coeff = sut.FitEval.model.constraints[sut.FitEval.model.TOTAL_PROTEIN_CONSTRAINT_ID].get_linear_coefficients([rxn.forward_variable])[rxn.forward_variable]
    assert sut.FitEval.model.constraints[sut.FitEval.model.TOTAL_PROTEIN_CONSTRAINT_ID].ub == tot_prot-intercept*1e3
    assert slope*1e3 == coeff


