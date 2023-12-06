"""
Minimal example for a fitness function class

Objective: Maximize the number of reaction deletions supporting a minimal growth rate

"""

import random

import deap.base
import numpy as np
import pandas as pd
from pathlib import Path
import os
from os.path import dirname, abspath

# set standard paths
FILE_PATH = Path(abspath(dirname(__file__)))
DATA_PATH = FILE_PATH.parents[0].joinpath("Data")

# seed random number generator
random.seed()

class FitnessEvaluation():
    #from kcat population (Bar-even, et al. 2011)
    KCAT_MU = 13.7
    KCAT_SIGMA = 1e2
    NUM_KCATS = 5



    def __init__(self, model=None, fixed_attr_list=[], processes=2, objective_id=str(),
                 valid_data_df = pd.DataFrame(), sigma_denominator:int = 10,
                 substrate_uptake_rates = [0.7, 11.3], substrate_uptake_id = 'EX_glc__D_e'):
        """Initialize fitness evaluation class for a genetic algorithm
        
        Inputs:
            :param cobra.core.Model model: Metabolic model in COBRA format
            :param list fixed_attr_list: Identifiers of attributes not to be used as solution variables
            :param int processes: Number of workers available (unused here)
            :param str objective_id: identifier of the objective function
            :param pd.DataFrame valid_data_df: Dataframe which contains the data to validate the results to
            :param int sigma_denominator: the factor determining the spread of the normal distribution from which
                    new kcat values will be sampled (kcat/sigma_denominator)
            :param list(float) substrate_uptake_rates: list with the substrate uptake rates to consider for
                    calculating the fitness (R^2 relative to the measurements at these substrate uptake rates)
            :param str substrate_uptake_id: identifier of the substrate uptake rate as defined in the model

        
        """
                    
        # save list of fixed attributes
        self.fixed_attr_list = fixed_attr_list

        #initiate reference values
        valid_data_df = valid_data_df.round({substrate_uptake_id: 5}) #round the substrate uptake ids for easy matching
        values_to_select = substrate_uptake_rates
        # check if the values are close, because of floating point errors this can deviate
        self.ref_data = valid_data_df[
            np.isclose(valid_data_df[substrate_uptake_id+'_ub'], values_to_select[0]) | np.isclose(
                valid_data_df[substrate_uptake_id+'_ub'],
                values_to_select[1])]


        # only get exchanges and growth rate
        self.growth_rate = [data for data in valid_data_df.columns if data.split('_')[0] == objective_id]
        self.reactions_with_data = [data for data in valid_data_df.columns if model.reactions.has_id(data)]

        # set the proper identifiers
        self.substrate_uptake_id = substrate_uptake_id
        self.substrate_uptake_rates = substrate_uptake_rates

        # set the factor determining the spread of the normal distribution from which new kcat values will be sampled
        self.sigma_denominator = sigma_denominator

        # initialize metabolic model
        self.init_model(model)
        
        
    ##########################################################################
    # NECESSARY CUSTOM FUNCTION
    ##########################################################################
    
    def init_model(self, model=None):
        """Pass and save the metabolic model
        Also constrain it to the condition for which we want to calculate the error between simulations and experiments
        Inputs:
            :param cobra.core.Model model: Metabolic model in COBRA format
        
        """
        self.model = model
        
    def compute_individuals_properties(self, pop) -> list:
        """Compute custom properties of a population's individuals. For each 
        individual return a dictionary with properties as keys
        
        Inputs:
            :param list pop: Individuals of a population
         
        Outputs:
            :param list ind_props: Property values for each individual. Return
                                    an empty list if no properties should be computed
        
        """
        ind_props = []
        # number of kcats changed
        for individual in pop:
            # determine which kcats have been changed
            kcats_changed = self._determine_changed_kcats(individual)
            # save property
            ind_props.append({
                'number_kcat_changes': len(kcats_changed)
                })
        
        return ind_props

    def _determine_changed_kcats(self, individual):
        changed = []
        for kcat, enzyme, reaction in zip(individual.kcat_list, individual.enzymes_to_eval, individual.reactions):
            kcat_init = self.model.enzymes.get_by_id(enzyme).rxn2kcat[reaction]['f']
            if kcat != kcat_init:
                changed.append(kcat)
        return changed
    
    ##########################################################################
    # DEAP-RELATED FUNCTIONS
    ##########################################################################
    def attribute_generator(self, probability=0, kcat_list = []) -> int:
        """Generates an attribute of a DEAP individual and is part of the 
        mutation operator
        
        Inputs:
            :param float probability: probability for returning a "0" as attr
        
        Outputs:
            :param int attr: representation of an attribute
        
        """

        if kcat_list == []:
            kcats = np.random.lognormal(mean = np.log10(self.KCAT_MU), sigma = np.log10(self.KCAT_SIGMA), size = self.NUM_KCATS)
            # individual.kcat_list = list(kcats)
            return kcats

        for i, kcat in enumerate(kcat_list):
            #use the sensitivity as mutation probability
            if random.random() <= probability:
                #scale needs to be positive, so prevent negative values
                # loc = mean, scale = sd, sd is defined as kcat/10 to make sure sd is in the same order of magntitude
                # as the kcat is
                new_kcat = self._mutate_kcat_value(kcat)
                kcat_list[i] = abs(new_kcat) #new kcat value should always be positive

        return kcat_list
    

    def init_fitness(self) -> dict:
        """Specifies parameters to set up a fitness function in DEAP
        - weights of fitness functions
        
        Outputs:
            :param dict param_fit: parameters to set up a DEAP fitness function
        
        """
        
        param_fit = {}
        # defines the weighting of the fitness functions
        # negative/positive weights minimize/maximize the objective function
        param_fit["weights"] = (1,) # must be a tuple so single and multi objective functions can be treated equally
        
        return param_fit
        
    
    def init_attribute(self, enzymes_to_eval:list) -> dict:
        """Specifies parameters to describe attributes of individuals in DEAP
        - number of attributes (returned)
        - type of an attribute (metabolic gene, regulator etc.)
        
        Each attribute entry (or solution variable) in 'self.individual_attr_list'
        is a dict specifying the identifier of the model variable (key: id) and
        the type of the variable (key: type). The order and numbers of attributes
        in the list matches an individual.
        
        Required format for "self.individual_attr_list":
            [{'id': 'var_1_id', 'type': 'var_1_type'}
             {'id': 'var_2_id', 'type': 'var_2_type'}
             ...
             ]
                
        
        Inputs:
            :param list enzymes_to_eval: list of enzyme ids which are to be mutated
            
        Outputs:
            :param dict param_attr: 
        
        """
        
        param_attr = {}
        
        # save attribute details for matching attributes to its model identity
        # here each model reaction is an attribute (solution variable)
        self.individual_attr_list = [
            {'id': enzyme, 'type': 'mutation'}
            for enzyme in enzymes_to_eval
            ]
        
        
        # number of attribute per individual
        param_attr["number_attributes"] = len(self.individual_attr_list)
        
        
        return param_attr
    
    
    
    
    def init_individual(self) -> dict:
        """Set up charactersitics and parameters of individuals in DEAP
        
        Outputs:
            :param dict param_ind: parameters for setting up an individual 

        """
        
        param_ind = {}
        

        # create individuals represnting metabolic genes as variables/targets
        # creator.create("Individual", list, fitness=creator.FitnessMin)
        param_ind["individual_type"] = list
        
        return param_ind
            
    
    
    ##########################################################################
    # FITNESS FUNCTION
    ##########################################################################
    def eval_fitness(self, individual) -> int:
        """Evaluate the fitness of an individual. It returns a tuple of 
        objective values.
        
        minimize the error to experimental measurements
        
        Minimize the computational workload in this function, since it will be
        frequently evaluated
        
        Inputs:
            :param list individual: a list of variables representing an individual (solution)
            :param PAModelpy.PAModel: model: copy of the model used in this iteration

        Outputs:
            :param int fitness: fitness of the current individual evaluated
        
        """

        # Here, the objective is to minimize the error to experimental measurements
        # (one low substrate uptake rate and a high substrate uptake rate are recommended)

        # apply kcat changes and compute metabolic functionalities
        fitness_list = []
        kcat_old = [self.model.enzymes.get_by_id(enz_id).get_kcat_values(individual.reactions[i]) for i, enz_id in enumerate(individual.enzymes_to_eval)]
        with individual.model as model:
            # change the kcat value of the enzymes with the highest sensitivity
            for i, enz_id in enumerate(individual.enzymes_to_eval):
                rxn = model.reactions.get_by_id(individual.reactions[i])
                model.change_kcat_value(enz_id,
                                  {individual.reactions[i]:
                                       {'f':individual.kcat_list[i]/(3600*1e-6), 'b':individual.kcat_list[i]/(3600*1e-6)}})
                model.constraints[f'EC_{enz_id}_f'].set_linear_coefficients({rxn.forward_variable:1/individual.kcat_list[i]}) #TODO remove after PAModelpy update
                model.constraints[f'EC_{enz_id}_b'].set_linear_coefficients(
                    {rxn.reverse_variable: 1 / individual.kcat_list[i]})
            # perform simulations
            for rate in self.substrate_uptake_rates:

                model.change_reaction_bounds(self.substrate_uptake_id,
                                                  lower_bound = 0, upper_bound = rate)
                model.slim_optimize(error_value=0)
                if model.solver.status != 'optimal':error = 1e2
                # calculate fitness (sum of simulation error to reactions with data)
                else: error = self._calculate_simulation_error(model, substrate_uptake=rate)

                fitness_list += [error]

            #average fitness:
            fitness = float(np.mean(fitness_list))
            individual.r_squared = fitness
            individual.fitness.values = [fitness]


        #revert kcat_changes
        for i, enz_id in enumerate(individual.enzymes_to_eval):
            model.change_kcat_value(enz_id,
                                     {individual.reactions[i]:
                                          {'f': kcat_old[i], 'b': kcat_old[i]}})
        # return a tuple of one element
        return tuple([float(fitness)])
    
    
    
    ##########################################################################
    # CUSTOM HELPER FUNCTIONS
    ##########################################################################
    
    def _calculate_simulation_error(self, model, substrate_uptake:float):
        error = []
        for rxn in self.reactions_with_data + self.growth_rate:
            #only select the rows which are filled with data

            ref_data_rxn = self.ref_data.dropna(axis = 0, subset = rxn)
            #if there are no reference data points, continue to the next reaction
            if len(ref_data_rxn) == 0: continue

            #select the right substrate uptake rate
            ref_data_rxn = ref_data_rxn[np.isclose(ref_data_rxn[self.substrate_uptake_id + '_ub'] ,substrate_uptake)|
                np.isclose(ref_data_rxn[self.substrate_uptake_id], substrate_uptake)]

            # calculate difference between simulations and validation data
            ref_data_rxn = ref_data_rxn.assign(simulation=model.reactions.get_by_id(rxn).flux)

            # error: squared difference
            ref_data_rxn = ref_data_rxn.assign(error=lambda x: abs((x[rxn] - x['simulation'])))
            error += [ref_data_rxn.error.mean()]

        return sum(error)/len(error)
        
    def _mutate_kcat_value(self, kcat:float, sensitivity:float = 0.5, toolbox:deap.base.Toolbox = None) -> float:
        """
        Helper function to mutate a kcat value based on a user-defined standard deviation scale.
        Can be used both with the deap toolbox (also need to include sensitivity to determine mutation probability),
        as well as without. In the latter case, a new kcat value will be sampled from np.random.normal Gaussian
        distribution.

        :param float kcat: kcat value to mutate
        :param float sensitivity: sensitivity of the enzyme related to the kcat to determine mutation probability
        :param deap.base.toolbox toolbox: deap toolbox

        :return: float new_kcat: mutated kcat value (sampled from normal distribution)
        """
        stdev = kcat/self.sigma_denominator
        if toolbox is not None:
            # mutate an individual with a mutation rate based on the sensitivity of the individual enzymes
            # the new value is samples from a gaussian distribution with mu being the original kcat value and
            # sigma being related to the kcat value to stay in sync with the order of magnitude of the original kcat
            return toolbox.mutate([kcat], mu=kcat, sigma=stdev, indpb=(1 - sensitivity))[0][0]
        else:
            # scale needs to be positive, so prevent negative values
            # loc = mean, scale = sd, sd is defined as kcat/10 to make sure sd is in the same order of magntitude
            # as the kcat is
            return float(np.random.normal(loc=kcat, scale=stdev))