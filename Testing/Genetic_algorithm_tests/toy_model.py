import sys
import os

from PAModelpy.Enzyme import Enzyme

import pandas as pd
sys.path.append('../../Modules/PAM_parametrization/genetic_algorithm_parametrization')
from core_parametrization import GAPO as GA

import cobra
from pathlib import Path

sys.path.append('../../Scripts/')
from pam_generation import setup_toy_pam

# Reference data from toy model simulations:
DATA_DIR = os.path.join(os.getcwd(), 'Data')
RESULT_DF_FILE = os.path.join(DATA_DIR, 'toy_model_simulations_ga.csv')

valid_data_df = pd.read_csv(RESULT_DF_FILE)
#start kcat: [1, 0.5, 1, 0.5 ,0.45, 1.5]
#aimed end kcat: [1, 0.5, 5, 0.1, 0.25, 1.5]

pamodel1 = setup_toy_pam()
pamodel2 = pamodel1.copy()
# print(pamodel.enzymes)

# %% initialize genetic algorithm
ga = GA(
    model=pamodel1, # Metabolic model,
    enzymes_to_eval = {'E3':{'reaction':'R3','kcat':1, 'sensitivity':0.5}, #should become 5
                       'E4':{'reaction':'R4','kcat':0.5, 'sensitivity':0.2},#should become 0.1
                       'E5':{'reaction':'R5','kcat':0.45, 'sensitivity':0.1}},#should become 0.25
    r_squared = 0.2,
    fitness_class="Fitfun_params", # filename (or module) of the fitness function class
    mutation_probability = 0.5, # probability with which an individual (solution) is mutated in a generation
    mutation_rate = 0.01, # probability with which an attribute (e.g. gene) of an individual is mutated
    population_size = 10, # number of individuals (solution) per population
    crossover_probability = 0.8, # probability with which two indivduals/offsprings are crossed over
    number_generations = 20, # number of consecutive generations per gene flow event
    number_gene_flow_events = 2, # number of gene flow events, i.e. merging of multiple
                                 # populations independently evolved on parallel workers
    init_attribute_probability=0.001, # probability with which attributes of initial individuals are mutated
    processes = 2, # number of parallel workers
    time_limit = 600, # time limit in seconds
    filename_save="toy_model_parametrization", # filename for saving results after every gene flow event
    folderpath_save=Path(r"Results"), # path for saving results
    overwrite_intermediate_results=True, # if true, saved intermediate results are overwritten
    valid_df = valid_data_df,
    objective_id = 'R7',
    substrate_uptake_rates = [0.001,0.091],
    substrate_uptake_id = 'R1'
    )

# %% start optimization
if __name__ == "__main__":
    # start genetic algorithm
    ga.start()
