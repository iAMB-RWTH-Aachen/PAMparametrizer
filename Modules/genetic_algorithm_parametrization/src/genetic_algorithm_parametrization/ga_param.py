"""
Genetic algorithm (GA) framework for optimization problems involving metabolic 
models in COBRA format
- DEAP (https://github.com/DEAP/deap) is used as a base GA
- GA can be generalized to any optimization problem incorporating M-models

"""

import random
from time import time, strftime
import numpy as np

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
            ind.fitness.values = fit
            
        return pop
        


    # main genetic algorithm iterations
    def main(self, pop, toolbox, start_time, fitfun, sensitivities, fitness_dict={}, pop_id="") -> (list, dict):
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
        
        # Begin the evolution
        while (g < self.number_generations) and ((time()-start_time) < self.time_limit):
            # A new generation
            g = g + 1
            print("({3}) Population {0}: Generation {1}/{2} --".format(pop_id, g, self.number_generations, print_time()))
            
            
            # get best individual of population for elitism
            elite = pop[0]
            for ind in pop:
                if ind.fitness._wsum() < elite.fitness._wsum():
                    # new best solution, swap
                    elite = ind
            #make clones of the elite individuals
            elite = list(map(toolbox.clone, [elite]))
    
            
            # Select the next generation individuals for breeding
            # offspring = toolbox.select(pop, int(population_size * selection_rate))
            offspring = toolbox.select(pop, population_size-1)
            # Clone the selected individuals
            offspring = list(map(toolbox.clone, offspring))
                      
    
            # Apply crossover and mutation on the offspring
            # the input individuals are being changed themselves
            for child1, child2 in zip(offspring[::2], offspring[1::2]):
    
                # cross the kcat_lists of two individuals with probability CXPB
                if random.random() < self.crossover_probability:
                    child1.kcat_list, child2.kcat_list =toolbox.mate(child1.kcat_list, child2.kcat_list)

        
                    # fitness values of the children
                    # must be recalculated later
                    del child1.fitness.values
                    del child2.fitness.values
    
    
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
                    # new_kcats =[toolbox.mutate([mutant.kcat_list[i]], mu=mutant.kcat_list[i], sigma = mutant.kcat_list[i]/10,
                    #                indpb=(1-sens))[0][0] for i, sens in enumerate(mutant.sensitivities)]
                    #make sure kcat is always positive
                    mutant.kcat_list= [kcat if kcat>0 else 0.001 for kcat in new_kcats]
                    for i in range(len(new_kcats)):
                        mutant[i] = mutant.kcat_list[i]
                    # mutant = [kcat if kcat>0 else 0.001 for kcat in new_kcats]
                    del mutant.fitness.values
        
            # Evaluate the individuals with an invalid fitness
            invalid_ind = []
            for ind in offspring:
                if not ind.fitness.valid:
                    # check if individual was already encountered
                    ind_hashable = ind_str(ind.kcat_list)
                    if ind_hashable in fitness_dict:
                        # load stored individual fitness
                        ind.fitness.values = fitness_dict[ind_hashable]
                    else:
                        # mark individual as invalid
                        invalid_ind.append(ind)

            
            fitnesses = map(toolbox.evaluate, invalid_ind)
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit
                # store fitness
                fitness_dict[ind_str(ind.kcat_list)] = fit

            
            # The population is entirely replaced by the offspring and the elite
            offspring = elite + offspring
            pop[:] = offspring
            
            # Gather all the fitnesses in one list and print the stats
            fits = [ind.fitness._wsum() for ind in pop]
            
            length = len(pop)
            mean = sum(fits) / length
            sum2 = sum(x*x for x in fits)
            std = abs(sum2 / length - mean**2)**0.5
            
            print("({4}) Population {0}: Generation {1}/{2} evaluated {3} individuals".format(
                pop_id, g, self.number_generations, len(invalid_ind), print_time())
                )
            print("\tMin %s" % min(fits), "\tMax %s" % max(fits), "\tAvg %s" % mean, "\tStd %s" % std)

        print("({1}) Population {0}: End of (successful) evolution --".format(pop_id, print_time()))
        
    
        
        return (pop, fitness_dict)
