"""
Genetic algorithm (GA) framework for optimization problems involving metabolic 
models in COBRA format
- DEAP (https://github.com/DEAP/deap) is used as a base GA
- GA can be generalized to any optimization problem incorporating M-models

"""

import random
from time import time, strftime

import deap.base
import numpy as np
from copy import deepcopy

# anonymus function for printing the time
print_time = lambda : strftime("%d/%m %H:%M:%S")



class Genetic_Algorithm():
    
    def __init__(self, crossover_probability=0.8, mutation_probability=0.5,
                 number_generations=20, time_limit=600):
        
        # number of generation run on a worker per gene flow event
        self.number_generations = number_generations
        # total time limit in seconds
        self.time_limit = time_limit
        # probability with which an individual is mutated
        self.mutation_probability = mutation_probability
        # probablity with which two offspring individuals are crossed over
        self.crossover_probability = crossover_probability
        
    
            
    def init_pop(self, toolbox, population_size,  evaluate_fitness=True) -> list:
        """Create an initial population of individuals (solutions)
        
        Inputs:
            :param deap.base.Toolbox toolbox: DEAP's toolbox class
            :param int population_size: Number of individuals
            :param bool evaluate_fitness: Evaluate the fitness of all individuals in the population
            
        Outputs:
            :param list pop: Individuals of a population
        
        
        """
        
        # create an initial population
        pop = toolbox.population(n=population_size)
        if evaluate_fitness:
            # evaluate the fitness of the population
            pop = self.evaluate_pop(pop, toolbox)
  
        return pop
    
    
    def evaluate_pop(self, pop, toolbox) -> list:
        """Evaluate fitness of all individuals within a population
        
        Inputs:
            :param list pop: Individuals of a population
            :param deap.base.Toolbox toolbox: DEAP's toolbox class
            
        Outputs:
            :param list pop: Individuals of a population
        """
        # Evaluate the entire population
        fitnesses = list(map(toolbox.evaluate, pop)) #[model for i in range(len(pop))]))
        for ind, fit in zip(pop, fitnesses):
            ind.fitness = fit
            # ind.fitness.values = fit
            
        return pop
        


    # main genetic algorithm iterations
    def main(self, pop, toolbox, start_time, fitfun, sensitivities, fitness_dict={}, pop_id="", print_progress = True) -> (list, dict):
        """Main genetic algorithm framework
        - supports elitism
        - mutation operator
        - selection and crossover operator
        - computing of fitness function values
        
        Inputs:
            :param list pop: Individuals of a population
            :param deap.base.Toolbox toolbox: DEAP's toolbox class
            :param time.time.time start_time: Start time of the genetic algorithm run
            :param FitnessEvaluation fitfun: fitness object to determine mutation distribution
            :param list sensitivities: list with the importance of each kcat towards changes in
                                        the growth rate in the upperlevel Protein Allocation Model
            :param dict fitness_dict: Fitness function values of individuals. Keys
                                        are hashable identifier of an individual
            :param str pop_id: Identifier of a population
            
        Outputs:
            :param list pop: Individuals of a population
            :param dict fitness_dict: Fitness function values of individuals. Keys
                                        are hashable identifier of an individual
        """
        
        # initialize random number generator
        random.seed()
          
        # initialize parameters
        population_size = len(pop)
        
        # Variable keeping track of the number of generations
        g = 0
        
        # lambda function to create a hashable representation of an individual
        ind_str = lambda ind: "".join([str(a) for a in ind])
        elite_fitness_prev = 0
        elite_kcat_prev =[0,0,0]
        
        # Begin the evolution
        while (g < self.number_generations) and ((time()-start_time) < self.time_limit):
            # A new generation
            g = g + 1

            if print_progress:
                print("({3}) Population {0}: Generation {1}/{2} --".format(pop_id, g, self.number_generations, print_time()))
            #
            # kcat_old = elite_kcat_prev.copy()
            # if g>1:
            #     print('elite in pop', elite[0] in pop, elite[0].fitness._wsum())
            #     for ind in pop:
            #         print(ind.fitness._wsum(), ind.kcat_list)
            # get best individual of population for elitism
            elite = self._get_best_individual_from_population(pop)
            # print('new elite fitness', elite.fitness._wsum())
            #make clones of the elite individuals
            elite = self._clone_elite([elite], toolbox)

            # import pytest
            # if elite[0].kcat_list == elite_kcat_prev:
            #     print('it is the same indiv!')
            #     if elite[0].fitness._wsum() == pytest.approx(elite_fitness_prev, abs = 1e-3):
            #         print('and the fitness is also the same!')
            # elif elite[0].fitness._wsum() < elite_fitness_prev:# & g>1:
            #     print('something fishy is happening here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            #
            #     print(elite[0].kcat_list == elite_kcat_prev, elite[0].kcat_list, elite_kcat_prev)
            elite_fitness_prev = elite[0].fitness._wsum()
            elite_kcat_prev = elite[0].kcat_list
            # print('elite fitness', elite[0].fitness._wsum())


            # Select the next generation individuals for breeding
            # offspring = toolbox.select(pop, int(population_size * selection_rate))
            # pop_elite = self._clone_elite(elite, toolbox)
            # if elite[0] in pop: pop = [item for item in pop if item!= elite[0]]#[item if item != elite[0] else pop_elite for item in pop]
            offspring = toolbox.select(pop, population_size-1)
            # Clone the selected individuals
            offspring = list(map(toolbox.clone, offspring))

            # Apply crossover and mutation on the offspring
            # the input individuals are being changed themselves
            for child1, child2 in zip(offspring[::2], offspring[1::2]):
    
                # cross the kcat_lists of two individuals with probability CXPB
                if random.random() < self.crossover_probability:
                    # child1.kcat_list, child2.kcat_list =toolbox.mate(child1.kcat_list, child2.kcat_list)
                    child1.kcat_list, child2.kcat_list =toolbox.mate(child1.kcat_list, child2.kcat_list)

                    # fitness values of the children
                    # must be recalculated later
                    del child1.fitness.values
                    del child2.fitness.values

            # print('after crossover', elite[0].fitness._wsum() == elite_fitness_prev,elite_fitness_prev, elite[0].fitness._wsum())
            # mutation operator
            for mutant in offspring:
                # if infeasible solutions are likely, mutate an individual only by chance
                if random.random() < self.mutation_probability:
                    # mutate an individual with a mutation rate based on the sensitivity of the individual enzymes
                    # the new value is samples from a gaussian distribution with mu being the original kcat value and
                    # sigma being related to the kcat value to stay in sync with the order of magnitude of the original kcat
                    new_kcats = [fitfun._mutate_kcat_value(kcat=kcat,
                                                                    sensitivity=sens,
                                                                    toolbox=toolbox)
                                 for kcat, sens in zip(mutant.kcat_list, sensitivities)]
                    #make sure kcat is always positive
                    mutant.kcat_list= [kcat if kcat>0 else 0.001 for kcat in new_kcats]
                    for i in range(len(new_kcats)):
                        mutant[i] = mutant.kcat_list[i]
                    del mutant.fitness.values
            # print('after mutation', elite[0].fitness._wsum() == elite_fitness_prev)
            invalid_ind = []
            for ind in offspring:
                if not ind.fitness.valid:
                    # check if individual was already encountered
                    ind_hashable = ind_str(ind.kcat_list)
                    if ind_hashable in fitness_dict:
                        # load stored individual fitness
                        ind.fitness = fitness_dict[ind_hashable]

                        # ind.fitness.values = fitness_dict[ind_hashable]
                    else:
                        # mark individual as invalid
                        invalid_ind.append(ind)

            # Evaluate the individuals with an invalid fitness
            fitnesses = map(toolbox.evaluate, invalid_ind)
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness = fit

                # ind.fitness.values = fit
                # store fitness
                fitness_dict[ind_str(ind.kcat_list)] = fit


            # The population is entirely replaced by the offspring and the elite
            offspring = elite + offspring
            pop[:] = offspring
            # Gather all the fitnesses in one list and print the stats
            fits = [ind.fitness._wsum() for ind in pop]
            # print('kcats did not change',kcat_old == elite[0].kcat_list)
            length = len(pop)
            mean = sum(fits) / length
            sum2 = sum(x*x for x in fits)
            std = abs(sum2 / length - mean**2)**0.5

            if print_progress:
                print("({4}) Population {0}: Generation {1}/{2} evaluated {3} individuals".format(
                    pop_id, g, self.number_generations, len(invalid_ind), print_time())
                    )
                print("\tMin %s" % min(fits), "\tMax %s" % max(fits), "\tAvg %s" % mean, "\tStd %s" % std)

        if print_progress:
            print("({1}) Population {0}: End of (successful) evolution --".format(pop_id, print_time()))
        return (pop, fitness_dict)

    def _get_best_individual_from_population(self, population):
        elite = population[0]
        for ind in population:
            if ind.fitness._wsum() > elite.fitness._wsum():
                # new best solution, swap
                elite = ind
        return elite

    def _clone_elite(self, pop_to_clone: list, toolbox: deap.base.Toolbox) -> list:
        new_population = []
        for indiv in pop_to_clone:
            kcat_list = deepcopy(indiv.kcat_list)
            new_indiv = toolbox.individual()
            new_indiv.kcat_list = kcat_list
            new_indiv.fitness = toolbox.evaluate(new_indiv)

            # new_indiv.fitness.values = toolbox.evaluate(new_indiv)
            new_population.append(new_indiv)

        return new_population