# Genetic Algorithm for Protein Allocation Model Parameter estimation
Framework for using a genetic algorithm to optimize the parametrization of Protein Allocation Models (PAMS).
For more information about the terminology and principles of genetic algorithms for metabolic optimization refer to this [paper](https://doi.org/10.3390/metabo8020033).  
The [DEAP toolbox](https://github.com/DEAP/deap) ([DEAP documentation](https://deap.readthedocs.io/en/master/index.html)) is used for basic genetic algorithm functionalities. 
This framework in based on the Genetic-algorithm-suite for optimization problems in metabolic models, with as initial 
usecase the computation of a minimal set of metabolic genes preserving certain core metabolic functionalities.

## Software Structure
The genetic algorithm in this module is build up out of 3 main parts:

1. `core_parametrizer`: the object which connects all other parts of the genetic algorithm.
2. `ga_params`: the actual genetic algorithm functionality
3. `Evaluation/Fitfun_params`: all the functions to handle a specific use cases of the genetic algorithm. It allows for a tailored fitness function, mutation function, individual format and storage of results


