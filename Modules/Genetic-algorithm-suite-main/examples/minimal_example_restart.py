"""
Restart the minimal example for using a genetic algorithm for finding solutions 
to problems involving metabolic models

- The optimization started here ties to a previous solution computed in "minimal_example.py"
- Here, the mutation rate is lowered to exploit the narrowed solution space
  from the previous run


RUN "minimal_example.py" BEFORE RUNNING THIS CODE FOR GENERATING INITIAL SOLUTIONS

"""

from genetic_algorithm_suite.core import GAMO as GA

import cobra
from pathlib import Path


# %% Load model
model = cobra.io.load_json_model('Models/iML1515.json')


# %% initialize genetic algorithm
ga = GA(
    model=model, # Metabolic model
    fitness_class="Fitfun_Minimal", # filename (or module) of the fitness function class
    mutation_probability = 0.5, # probability with which an individual (solution) is mutated in a generation
    mutation_rate = 0.001, # probability with which an attribute (e.g. gene) of an individual is mutated
    population_size = 10, # number of individuals (solution) per population
    crossover_probability = 0.8, # probability with which two indivduals/offsprings are crossed over 
    number_generations = 20, # number of consecutive generations per gene flow event
    number_gene_flow_events = 2, # number of gene flow events, i.e. merging of multiple 
                                 # populations independently evolved on parallel workers
    init_attribute_probability=0.001, # probability with which attributes of initial individuals are mutated 
    processes = 2, # number of parallel workers
    time_limit = 600, # time limit in seconds
    filename_save="minimal_example", # filename for saving results after every gene flow event
    folderpath_save=Path(r"Results"), # path for saving results
    overwrite_intermediate_results=True # if true, saved intermediate results are overwritten
    )


# %% restart optimization from a previous, intermediate results
if __name__ == "__main__":
    previous_run_path = Path("Results/minimal_example.json")
    ga.restart(previous_run_path)