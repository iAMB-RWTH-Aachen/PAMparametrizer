"""
Start computation of a minimal set of metabolic genes, which preserves a custom
set of metabolic functionalities. For the protected functionalities refer to
"Metabolic_Genes_Reduction.py" and the "_compute_metabolic_functionalities" function

THE FITNESS FUNCTION IS COMPUTATIONALLY VERY INTENSIVE. IT IS RECOMMENDED TO
RUN THIS CODE FOR >24h ON MULTIPLE WORKERS  (>20) AND RESTART THE OPTIMIZATION
WITH A DECREASED MUTATION RATE.

For re-starting the optimization from a previous population/solution with 
altered parameters refer to the minimal example

To set up custom fitness function classes refer to "Fitfun_Minimal.py" (in
"\src\genetic_algorithm_suite\Evaluation") for a template and minimum requirements.

Custom fitness function classes may be stored in "\src\genetic_algorithm_suite\Evaluation"
and used by parsing the filename via the "fitness_class" keyword argument (see below),
or the respective module is loaded and parsed via the same keyword argument.

"""

from genetic_algorithm_suite.core import GAMO as GA

import cobra
from pathlib import Path


# %% Load model
model = cobra.io.load_json_model('Models/iML1515.json')


# %% initialize genetic algorithm
ga = GA(
    model=model, # Metabolic model
    fitness_class="Metabolic_Genes_Reduction", # filename (or module) of the fitness function class
    mutation_probability = 0.5, # probability with which an individual (solution) is mutated in a generation
    mutation_rate = 0.01, # probability with which an attribute (e.g. gene) of an individual is mutated
    population_size = 10, # number of individuals (solution) per population
    crossover_probability = 0.8, # probability with which two indivduals/offsprings are crossed over 
    number_generations = 20, # number of consecutive generations per gene flow event
    number_gene_flow_events = 2, # number of gene flow events, i.e. merging of multiple 
                                 # populations independently evolved on parallel workers
    init_attribute_probability=0.001, # probability with which attributes of initial individuals are mutated 
    processes = 2, # number of parallel workers
    time_limit = 3600, # time limit in seconds
    filename_save="ga_genome_reduction", # filename for saving results after every gene flow event
    folderpath_save=Path(r"Results"), # path for saving results
    overwrite_intermediate_results=False # if true, saved intermediate results are overwritten
    )

# %% start optimization
if __name__ == "__main__":   
    # start genetic algorithm
    ga.start()

