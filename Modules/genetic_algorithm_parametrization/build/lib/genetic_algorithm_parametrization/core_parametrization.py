"""
Genetic algorithm (GA) for the prediction genome reduction paths based on metabolic models
- DEAP (https://github.com/DEAP/deap) is used as a base GA
- GA can be generalized to any optimization problem incorporating M-models

"""

# Disable gurobi logging output
try:
    import gurobipy
    gurobipy.setParam("OutputFlag", 0)
except ImportError:
    pass

import random
from time import time
from time import strftime
import pickle
from pathlib import Path
import importlib
import pandas as pd
import json

from multiprocessing import Pool

from deap import base
from deap import creator
from deap import tools

from genetic_algorithm_parametrization.ga_param import Genetic_Algorithm

# anonymus function for printing the time
print_time = lambda : strftime("%d/%m %H:%M:%S")


class GAPO():
    
    def __init__(self, model=None,
                 enzymes_to_eval: dict = {},  #dict of enz.id:{reaction, kcat, sensitivity}
                 r_squared: float = 1,
                 fitness_class = "Fitfun_params",
                 mutation_probability=0.5, mutation_rate=0.05, population_size=30,
                 crossover_probability=0.8, number_generations=20, number_gene_flow_events=10,
                 processes=2, time_limit=600, init_attribute_probability=0,
                 fixed_attributes=[],
                 folderpath_save=Path("Results"), filename_save="ga_results",
                 overwrite_intermediate_results=True,
                 objective_id = 'BIOMASS', valid_df = pd.DataFrame(),
                 sigma_denominator:int=10,
                 substrate_uptake_rates = [0.7,11.3], substrate_uptake_id = 'EX_glc__D_e'):
        
        if not model:
            self.model = model
        else:
            self.model = model.copy()

        #store information for initialization of the population
        self.rxns = list()
        self.kcat_list = list()
        self.sensitivity_list = list()
        self.enzymes_to_eval = list()
        for enzyme_id, values in enzymes_to_eval.items():
            self.enzymes_to_eval += [enzyme_id]
            self.rxns += [values['reaction']]
            self.kcat_list += [values['kcat']]
            self.sensitivity_list += [values['sensitivity']]


        # Specify GA parameters
        # probability for mutating an individual
        self.mutation_probability = mutation_probability
        # probability for mutating an attribute
        self.mutation_rate = mutation_rate
        # number of individuals per population
        self.population_size = population_size
        # probability for mutating an individual
        self.crossover_probability = crossover_probability
        # number of generations
        self.number_generations = number_generations
        # number of gene drifts
        self.number_gene_flow_events = number_gene_flow_events
        # probability of generating a "0" as attribute for initial population
        self.init_attribute_probability = init_attribute_probability
        # number of parallel processes
        self.processes = processes
        # maximum time limit [sec]
        self.time_limit = time_limit
        # filename for saving results
        self.filename_save = filename_save
        # intermediate results are overwritten if true
        # if false, results after a gene flow event are saved separately 
        self.overwrite_intermediate_results = overwrite_intermediate_results
   
        
        # summarize parameters in dictionary
        self.ga_parameters = {
            'current_gene_flow_number': 0,
            'mutation_probability': mutation_probability,
            'mutation_rate': mutation_rate,
            'population_size': population_size,
            'crossover_probability': crossover_probability,
            'number_generations': number_generations,
            'number_gene_flow_events': number_gene_flow_events,
            'init_attribute_probability': init_attribute_probability,
            'processes': processes,
            'time_limit': time_limit,
            'filename_save': filename_save,
            'folderpath_save': str(folderpath_save),
            'overwrite_intermediate_results': overwrite_intermediate_results,
            }
        

        # setup save folder
        self.folderpath_save = folderpath_save
        if not folderpath_save.is_dir():
            folderpath_save.mkdir()


        # load genetic algorithm
        self.ga = Genetic_Algorithm(
            crossover_probability = self.crossover_probability,
            mutation_probability = self.mutation_probability,
            number_generations = self.number_generations,
            time_limit = self.time_limit
            )


        # load preinstalled or use parsed custom fitness function evaluation class
        if isinstance(fitness_class, str):
            # load preinstalled module
            self.fitness_class = importlib.import_module("genetic_algorithm_parametrization.Evaluation."+fitness_class)
        else:
            self.fitness_class = fitness_class



        # Set up fitness evaluation class
        self.FitEval = self.fitness_class.FitnessEvaluation(
            model=self.model,
            processes=processes,
            fixed_attr_list=fixed_attributes,
            valid_data_df = valid_df,
            sigma_denominator = sigma_denominator,
            objective_id = objective_id,
            substrate_uptake_rates = substrate_uptake_rates,
            substrate_uptake_id = substrate_uptake_id)

        
        self._init_deap_fitness() # initialize the fitness function
        self._init_deap_individual(r_squared) # initialize deap individual representation


    # execute genetic algorithm
    def start(self):

        # initialize timing
        start_time = time()
    
    
        # initialize DEAP toolbox
        print("({}) Initialize DEAP toolbox --".format(print_time()))

        self.toolbox = self._init_deap_toolbox() # initialize the toolbox
        
        # save evaluation class
        with open(self.folderpath_save.joinpath(self.filename_save+".pickle"), "wb") as f:
            pickle.dump(self.FitEval, f)
       
        # initialize populations (multiprocessing)
        print("({}) Initialize population --".format(print_time()))
            
        with Pool(processes=self.processes) as pool:
            pops = pool.starmap(
                self.ga.init_pop,
                [(self.toolbox,self.population_size,True)
                 for i in range(self.processes)]
                )

        # start optimization with parallel gene flow events
        print("({}) Start optimization --".format(print_time()))
        self.pops_final = self._parallel_gene_flow(pops, self.toolbox, start_time)

        # evaluate final population
        # results = self.FitEval.eval_population(sum(self.pops_final, []), self.filename_save)
        print("({}) Evaluate final population --".format(print_time()))
        self._save_population(sum(self.pops_final, []))
        
        
    
    def restart(self, filepath_previous_pop):
        """Restart genetic algorithm with the final population from a preceeding run
        
        Inputs:
            :param pathlib.Path filepath_previous_pop: Path to previous genetic algorithm results
            
        Outputs:
            
        
        
        """
        
        # helper functions
        def get_fitness(ind):
            return ind["fitness_weighted_sum"]
        
        
        # initialize timing
        start_time = time()
        
        print("({}) Load previous population data --".format(print_time()))
        # load previous final population data
        if isinstance(filepath_previous_pop, str):
            # get Path to file
            filepath_previous_pop = Path(filepath_previous_pop)
            
        with open(filepath_previous_pop, "r") as f:
            pop_previous_dict = json.load(f)
            
        self.pop_previous_dict = pop_previous_dict 
        
        # initialize DEAP toolbox
        print("({}) Initialize DEAP toolbox --".format(print_time()))
        self.toolbox = self._init_deap_toolbox()
        
        # save evaluation class
        with open(self.folderpath_save.joinpath(self.filename_save+".pickle"), "wb") as f:
            pickle.dump(self.FitEval, f)
            
        # initialize population
        print("({}) Initialize population and populate with previous individuals --".format(print_time()))
        pops = [self.ga.init_pop(self.toolbox, self.population_size, evaluate_fitness=False) for p in range(self.processes)]

        # sort previous population (maximization of weighted fitness sum)
        pop_previous_dict["population"].sort(reverse=True, key=get_fitness)
        
        # extract best individuals from previous population
        if (self.population_size*self.processes) >= len(pop_previous_dict["population"]):
            pop_previous = pop_previous_dict["population"]
        else:
            pop_previous = pop_previous_dict["population"][:self.population_size*self.processes]
         
      

        # map previous attributes list to the current list
        mapping_vector = [None for a in self.FitEval.individual_attr_list]
        for j in range(len(self.FitEval.individual_attr_list)):
            attr_curr = self.FitEval.individual_attr_list[j]            
            for i in range(len(pop_previous_dict["attributes"])):
                if (pop_previous_dict["attributes"][i]["type"] == attr_curr["type"]) \
                    and (pop_previous_dict["attributes"][i]["id"] == attr_curr["id"]):
                    mapping_vector[j] = i
                    break
         

        # merge previous population into current population
        random.shuffle(pop_previous) # shuffle previous population
        
        
        i = 0 # counter individuals in previous population
        j = 0 # counter individuals in current sub-population
        while i < len(pop_previous):
            for p in range(len(pops)):
                # if there are no more individuals from the previous run left to allocate, break
                if i >= len(pop_previous):
                    break
                # match attributes of the previous with the current run
                for a in range(len(pops[p][j])):
                    if not(not(mapping_vector[a])):
                        pops[p][j][a] = pop_previous[i]["attribute_values"][mapping_vector[a]]
                        
                i += 1 
            j += 1
                
        
        # evaluate fitness of population
        print("({}) Evaluate fitness of population --".format(print_time()))
        with Pool(processes=self.processes) as pool:  
            pops = pool.starmap(self.ga.evaluate_pop, [(pops[i], self.toolbox) for i in range(self.processes)])

        # re-start optimization
        # update number of gene flow events
        self.number_gene_flow_events = self.number_gene_flow_events + pop_previous_dict['ga_parameters']['current_gene_flow_number']
        self.ga_parameters['number_gene_flow_events'] = self.number_gene_flow_events
        
        print("({}) Start optimization --".format(print_time()))
        self.pops_final = self._parallel_gene_flow(
            pops,
            self.toolbox,
            start_time,
            previous_drifts=pop_previous_dict['ga_parameters']['current_gene_flow_number'])

        # evaluate and save final population
        print("({}) Evaluate final population --".format(print_time()))
        self._save_population(sum(self.pops_final, []))
        


    def _init_deap_toolbox(self):
        
        #  initialize DEAP toolbox
        #           this is the central module incorporating all the information, param
        #           parameters, and data of the genetic algorithm problem
        toolbox = base.Toolbox()
        
        # individual generator: mutating the list of kcat values
        toolbox.register("individual_generator", self.FitEval.attribute_generator,
                         self.init_attribute_probability)
        

        # register individual representation
        param_attr = self.FitEval.init_attribute(self.enzymes_to_eval)

        toolbox.register("individual", toolbox.individual_generator, creator.Individual)
        
        # define the population to be a list of individuals
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)
        
        #  Operator registration
        # register the goal / fitness function
        toolbox.register("evaluate", self.FitEval.eval_fitness)
        
        # register the crossover operator
        toolbox.register("mate", tools.cxTwoPoint)
        
        # register a mutation operator with a probability to
        # flip an attribute
        toolbox.register("mutate", tools.mutGaussian, indpb=self.mutation_rate)
        
        # operator for selecting individuals for breeding the next
        # generation: each individual of the current generation
        # is replaced by the 'fittest' (best) of three individuals
        # drawn randomly from the current generation.
        toolbox.register("select", tools.selTournament, tournsize=2)
    
        
        
        return toolbox
    
    
    def _init_deap_individual(self, r_squared):
        
        # create individuals representing metabolic genes as variables/targets
        param_ind = self.FitEval.init_individual()
        creator.create("Individual", param_ind["individual_type"],
                       fitness=creator.FitnessObj(), reactions = self.rxns, enzymes_to_eval = self.enzymes_to_eval,
                       kcat_list = self.kcat_list, sensitivities = self.sensitivity_list, r_squared = r_squared)
        
    def _init_deap_fitness(self):
        
        param_fit = self.FitEval.init_fitness()
        creator.create("FitnessObj", MyFitness, weights=param_fit["weights"])
    
    
    def _parallel_gene_flow(self, pops, toolbox, start_time, previous_drifts=0):
        
        # initialize fitness dictionary
        fitness_dict = {}
    
        drift = previous_drifts
        
        print('\nTime left:', '{0}min'.format(round((self.time_limit-time()+start_time)/60, 1)))
        
        while (drift < self.number_gene_flow_events) and ((time()-start_time) < self.time_limit):
            
            drift += 1
            print("\n-- Gene drift %i --" % drift)
            print('Time left:', '{0}min'.format(round((self.time_limit-time()+start_time)/60, 1)))  
            print("Create populations --")
            # shuffle individuals among populations
            unq_pop = sum(pops, []) # unique set of all individuals
            shuffled_idx = random.sample([i for i in range(len(unq_pop))], len(unq_pop))
            unq_ind_idx = 0
            for pop_idx in range(len(pops)):
                for ind_idx in range(len(pops[pop_idx])):
                    # randomly pick an individual from the total set of individuals
                    pops[pop_idx][ind_idx] = unq_pop[shuffled_idx[unq_ind_idx]]
                    unq_ind_idx += 1
                    
               
            # multiprocessing
            print("Start genetic algorithm --")
            with Pool(processes=self.processes) as pool:
                # distribute populations to separate workers
                gen_results = pool.starmap(self.ga.main, [(pops[i], toolbox, start_time, self.FitEval, fitness_dict, str(i+1)) for i in range(len(pops))])

                # extract populations
                print("Postprocess evolved population --")
                for i in range(len(gen_results)):
                    pops[i] = gen_results[i][0] # individuals of the returned population

                
            # evaluate and save all populations
            self.ga_parameters['current_gene_flow_number'] = drift # update current number of gene flow event
            if self.overwrite_intermediate_results:
                suffix = ''
            else:
                suffix = '_{}'.format(drift)
                
            self._save_population(sum(pops, []), suffix=suffix)
            
            print('Time left:', '{0}min'.format(round((self.time_limit-time()+start_time)/60, 1)))

        
        
        return pops
    
    
    def _save_population(self, pop, suffix=""):
        """Saves a population and the metadat for the genetic algorithm in
        Excel and machine readable .json format
        
        Inputs:
            :param list pop: Individuals of a population
            :param str suffix: Suffix for the filename
            
        Outputs:
        
        
        """
         
        # determine save path
        save_path = str(self.folderpath_save.joinpath(self.filename_save+suffix))
        
        # get attributes list from custom fitness evaluation class
        individual_attr_list = self.FitEval.individual_attr_list
        # gte fixed attributes list from custom fitness evaluation class
        fixed_attr_list = self.FitEval.fixed_attr_list
        
        # get best individual and save individual's parameters
        # get custom properties of a population's individuals
        pop_custom_properties = self.FitEval.compute_individuals_properties(pop)
        if len(pop_custom_properties) != len(pop):
            # property list doesn't match population
            pop_custom_properties = [{} for i in range(len(pop))]
            print('Custom population properties are invalid')
        
        pop_list = []
        best_ind = pop[0] # initialize best individual
        for ind, ind_custom_properties in zip(pop, pop_custom_properties):
            # get properties
            ind_properties = {
                "fitness_weighted_sum": ind.fitness._wsum(), # fitness value
                "attributes": ",".join([str(j) for j in ind.kcat_list]) # attributes list as string
                }
            # merge with custom properties and save
            pop_list.append({**ind_properties, **ind_custom_properties})
            
            # save best individual
            if best_ind.fitness._wsum() < ind.fitness._wsum():
                best_ind = ind
              
        pop_frame = pd.DataFrame(pop_list)    
        
        # sort population
        pop_frame = pop_frame.sort_values(by="fitness_weighted_sum", axis=0, ascending=False)
        
        # evaluate and save best individual
        best_ind_frame = pd.DataFrame({
            "id": [attr['id'] for attr in individual_attr_list],
            "type": [attr['type'] for attr in individual_attr_list],
            "value": [v for v in best_ind.kcat_list]
            })
        # sort frame
        best_ind_frame = best_ind_frame.sort_values(by="id", axis=0, ascending=True)
        
        
        # write attributes list
        attributes_frame = pd.DataFrame({
                "id": [attr['id'] for attr in individual_attr_list],
                "type": [attr['type'] for attr in individual_attr_list],
                })
        
        
        fixed_attributes_frame = pd.DataFrame({
                "id": fixed_attr_list,
                })
                
        # save metadata (genetic algorithm parameters)
        metadata_frame = pd.DataFrame({
            'parameter': self.ga_parameters.keys(),
            'value': self.ga_parameters.values(),
            })
                
                
        # save population and attributes list in machine readable form
        pop_dict = {"population": [],
                    "attributes": individual_attr_list, 
                    "fixed_attr_list": fixed_attr_list,
                    "fitness": {"weights": ind.fitness.weights},
                    "ga_parameters": self.ga_parameters
                    }
        
        for ind in pop:
            pop_dict["population"].append({
                "attribute_values": [a for a in ind.kcat_list],
                "fitness_weighted_sum": ind.fitness._wsum(),
                "fitness_values": ind.fitness.values
                })

        # save population
        # save frames (human-readable)
        writer = pd.ExcelWriter(save_path+'.xlsx', engine="openpyxl")   
        
        # metadata
        metadata_frame.to_excel(writer, sheet_name="metadata", index=False)
        # best individual
        best_ind_frame.to_excel(writer, sheet_name="best_individual", index=False)
        # whole population
        pop_frame.to_excel(writer, sheet_name="final_population", index=False)
        # attribute types and names
        attributes_frame.to_excel(writer, sheet_name="attributes_list", index=False)
        # fixed attribute types, values, and names
        fixed_attributes_frame.to_excel(writer, sheet_name="fixed_attributes_list", index=False)
        
        writer.save()
        
        # save dict (machine-readable)
        with open(save_path+'.json', "w") as f:
            json.dump(pop_dict, f)
        
    

    
# %% Override DEAP's base.Fitness class with more meaningful comparator methods
# DEAP's comparators work with Python's "funny" tuple comparator
 
class MyFitness(base.Fitness):
    
    def _wsum(self):
        # return sum of weighted values
        return sum(self.wvalues)
       
    def __le__(self, other):
        return self._wsum() <= other._wsum()

    def __lt__(self, other):
        return self._wsum() < other._wsum()

    def __eq__(self, other):
        return self._wsum() == other._wsum()           
        
            
                

    