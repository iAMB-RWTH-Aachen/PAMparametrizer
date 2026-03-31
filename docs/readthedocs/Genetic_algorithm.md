# Genetic Algorithm for Protein Allocation Model Parameter Estimation
This subdirectory contains the framework for using a genetic algorithm to optimize the parametrization of Protein Allocation Models (PAMs).
For more information about the terminology and principles of genetic algorithms for metabolic optimization refer to this [paper](https://doi.org/10.3390/metabo8020033).  
The [DEAP toolbox](https://github.com/DEAP/deap) ([DEAP documentation](https://deap.readthedocs.io/en/master/index.html)) is used for basic genetic algorithm functionalities. 
This framework in based on the Genetic-algorithm-suite for optimization problems in metabolic models, with as initial 
usecase the computation of a minimal set of metabolic genes preserving certain core metabolic functionalities.

## Software Structure
The genetic algorithm in this module is build up out of 3 main parts:

1. `core_parametrizer`: the object which connects all other parts of the genetic algorithm.
2. `ga_params`: the actual genetic algorithm functionality
3. `Evaluation/Fitfun_params`: all the functions to handle a specific use cases of the genetic algorithm. It allows for a tailored fitness function, mutation function, individual format and storage of results

When the `core_parametrizer` is called, it will initialize both a genetic algorithm instance and the fitness function.
The core object will initialize the DEAP toolbox and configure the individual. When the `start` function is called, 
the core will initialize the genetic algorithm with a population, optionally in a multithreaded manner. The core will
also start the genetic algorithm function in the genetic algorithm object and parse the results.

The `ga_params` file contains the `Genetic_Algorithm` (genetic algorithm for parameter optimization) object. The `GAPO`'s `main` 
function is called by the core. In the main function all the steps which are involved in the execution of the genetic 
algorithm are listed. This includes selection of the elite, crossover, mutation and evaluation of the individuals in a 
population.

In `Evaluation/Fitfun_params` the manually `FitnessEvaluation` class can be found. These class contains all the functions
which are required to 'personalize' the genetic algorithm for a specific usecase. Here there are functions for initialization
of several attributes, parsing of the properties which needs to be written to a file, mutating a kcat value, and 
importantly, to calculate the fitness of an individual.

## Genetic Algorithm workflow
The `main` function in the `Genetic_Algorithm` object follows the following steps (in pseudocode)

>`while (g < self.number_generations) and ((time()-start_time) < self.time_limit):`
>>1. **Elitism**: 
>>> * Select the best individual in a population 
>>> * Make a copy of the elite to retain the best solution
>>2. **Selection**: 
>>> * Randomly draw 3 individuals from the population
>>> * Select best individual of the 3 individuals
>>> * Replace 2 worse individual by copies of the best individual
>>3. **Crossover**:
>>> * Cross the lists of kcat values of two individuals with predefined crossover probability
>>4. **Mutation**:
>>> * Randomly draw a new kcat value from a distribution based on a mutation probability
>>> * If the model is likely to become infeasible, consider to decrease the mutation probability
>>4. **Fitness calculation**:
>>> * Calculate fitness for each individual which underwent some change
>>> * Replace old population with new population


## Fitness Evaluation adjustment
The path to a `FitnessEvaluation` object is one of the inputs for running the genetic algorithm. This allows the user to 
put in a personalized fitness function, kcat mutation method and to change the attributes which are saved. This 
repository already contains two: one FitnessEvaluation which mutates the kcat values using gaussian sampling given a 
specific standard deviation and one which uses uniform sampling, sampling from 0 to 2*kcat, whith the diffusion limit (1e6) providing an
upper bound for the maximal k<sub>cat</sub> value. Please note that the genetic algorithm modifies the k<sub>cat<\sub>
by changing the corresponding coefficient in the constraint relating the enzyme concentration to the reaction flux, without using
the utilities from PAModelpy. The kcat values are thus saved in the genetic algorithm as model coefficients (1/k<sub>cat</sub>, with k<sub>cat</sub> in 1/h).

By default the genetic algorithm uses uniform sampling, as this has been proven to yield the best results.

