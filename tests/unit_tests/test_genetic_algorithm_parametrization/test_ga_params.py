import unittest
import pytest
import os
import pandas as pd
from io import StringIO

from Scripts.Testing.Genetic_algorithm_tests.toy_model import init_toy_parametrization_ga
from Scripts.toy_ec_pam import evaluate_toy_model_fitness
from Scripts.pam_generation import setup_toy_pam




class TestGaParam(unittest.TestCase):
    def test_evaluate_pop_correct_fitness(self):
        """
        Evaluates the evaluate_pop function, but indirectly the toolbox.evaluate function which is then called
        :return:
        """
        # Arrange
        # Reference data from toy model simulations:
        DATA_DIR = os.path.join(os.getcwd(), 'Scripts', 'Testing', 'Data')
        RESULT_DF_FILE = os.path.join(DATA_DIR, 'toy_model_simulations_ga.csv')
        valid_data_df = pd.read_csv(RESULT_DF_FILE)
        #setup genetic algorithm and population
        toy_model_gaprot = init_toy_parametrization_ga(valid_data_df)
        toolbox = toy_model_gaprot._init_deap_toolbox()
        toy_ga = toy_model_gaprot.ga
        population = toolbox.population(n=1)

        # Act
        population = toy_ga.evaluate_pop(population, toolbox)
        fitness_simulated = population[0].fitness.values[0]

        # Assert
        #adjust for altered kcat_values
        toy_pam = setup_toy_pam(kcat_fwd = [1, 0.5]+population[0].kcat_list+[1.5])
        fitness_validation = evaluate_toy_model_fitness(toy_pam, reference_data_file_path = RESULT_DF_FILE)

        #1e-6 is solver feasibility tolerance
        assert fitness_validation == pytest.approx(fitness_simulated, 1e-6)


if __name__ == '__main__':
    unittest.main()
