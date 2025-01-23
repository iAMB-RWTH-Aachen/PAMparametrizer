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
from scipy.stats import linregress
from typing import Union, Literal
from ..core_parametrization_gaussian import MyFitness

from Modules.utils.error_calculation import calculate_r_squared_for_reaction
from Modules.utils.sector_config_functions import change_translational_sector_with_config_dict

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

    def __init__(self, model=None,
                 translational_sector_config:dict = None,
                 fixed_attr_list=[],
                 objective_id=str(),
                 valid_data = dict(),
                 sigma_denominator:int = 10,
                 error_weights: dict = {},
                 substrate_uptake_rates = {'EX_glc__D_e':[0.7, 11.3]},
                 substrate_uptake_id = 'EX_glc__D_e',
                 enzymes_to_evaluate: dict[Literal['reaction', 'kcat', 'sensitivity'], Union[str, dict, float]] = None):
        """Initialize fitness evaluation class for a genetic algorithm
        
        Args:
            cobra.core.Model: Metabolic model in COBRA format
            translational_sector_config (dict of dict): Dictionary with the slope and intercept of the translational
                sector configuration for each substrate.
                Format:
                {'substrate_uptake_id':{
                    'slope':float, #slope in g/mmol/h
                    'intercept':float #intercept in g/mmol
                    }
                }
            fixed_attr_list: Identifiers of attributes not to be used as solution variables
            processes: Number of workers available (unused here)
            objective_id: identifier of the objective function
            valid_data: dictionary of substrate uptake reaction and Dataframe which contains
                the data to validate the results to for that specific substrate
            sigma_denominator: the factor determining the spread of the normal distribution from which
                    new kcat values will be sampled (kcat/sigma_denominator)
            substrate_uptake_rates: dict with substrate id, list with the substrate uptake rates to consider for
                    calculating the fitness (R^2 relative to the measurements at these substrate uptake rates)
            substrate_uptake_id: identifier of the substrate uptake rate as defined in the model
            enzymes_to_evaluate: dictionary with enzyme-reaction-direction mapping as given to the genetic algorithm
                    enz.id:{reaction, kcat, sensitivity}

        
        """
                    
        # save list of fixed attributes
        self.fixed_attr_list = fixed_attr_list

        self.valid_data = valid_data
        self.growth_rate = {}
        self.reactions_with_data = {}
        self.substrate_uptake_rates = {}
        self.translational_sector_config = translational_sector_config

        # only get exchanges and growth rate
        for substr_uptake, valid_data_df in valid_data.items():
            self.growth_rate[substr_uptake] = [data for data in valid_data_df.columns if data.split('_')[0] == objective_id]
            self.reactions_with_data[substr_uptake] = [data for data in valid_data_df.columns if model.reactions.has_id(data)]
            self.substrate_uptake_rates[substr_uptake] = [round(rate, 6) for rate in substrate_uptake_rates[substr_uptake]]

        self.weights = error_weights
        self.enzymes_to_evaluate = enzymes_to_evaluate

        # set the proper identifiers
        self.substrate_uptake_id = substrate_uptake_id

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
        self.totprot_start_ub = self.model.constraints[self.model.TOTAL_PROTEIN_CONSTRAINT_ID].ub
        self.model.constraints[self.model.TOTAL_PROTEIN_CONSTRAINT_ID].lb = 0
        self.tps_intercept = self.model.sectors.get_by_id('TranslationalProteinSector').intercept


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
        for kcat, enzyme, direction, reaction in zip(individual.kcat_list, individual.enzymes_to_eval, individual.directions,individual.reactions):
            kcat_init = self.model.enzymes.get_by_id(enzyme).rxn2kcat[reaction][direction]
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
        
    
    def init_attribute(self, enzymes_to_eval:list, direction:list, rxns_to_eval:list) -> dict:
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
            :param list rxns_to_eval: list of reactions ids which the enzymes are catalyzing
            
        Outputs:
            :param dict param_attr: 
        
        """
        
        param_attr = {}
        
        # save attribute details for matching attributes to its model identity
        # here each model reaction is an attribute (solution variable)
        self.individual_attr_list = [
            {'id': enzyme, 'type': 'mutation', 'rxn_id': rxn, 'direction': dir}
            for enzyme, dir, rxn in zip(enzymes_to_eval, direction, rxns_to_eval)
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
        """Evaluate the fitness of an individual. It returns newly calculated fitnesses

        The fitness is defined as the R^2 value when comparing to experimental data

        
        Inputs:
            :param list individual: a list of variables representing an individual (solution)
        Outputs:
            :param int fitness: fitness of the current individual evaluated
        
        """
        self._reset_translational_sector()

        # apply kcat changes and compute metabolic functionalities
        self._change_kcat_values_for_individual(individual)

        # perform simulations and save results
        fluxes = {substr_uptake: pd.DataFrame(
            columns = ['substrate'] + self.reactions_with_data[substr_uptake]
        ) for substr_uptake in self.reactions_with_data.keys()}
        error = []

        for substrate_uptake_id, fluxes_df in fluxes.items():
            if self.translational_sector_config is not None:
                change_translational_sector_with_config_dict(self.model,
                                                             self.translational_sector_config[substrate_uptake_id],
                                                             substrate_uptake_id)
            for rate in self.substrate_uptake_rates[substrate_uptake_id]:
                new_row = [rate] + [0] * len(fluxes_df.columns[1:])
                fluxes_df.loc[len(fluxes_df)] = new_row
                if rate >= 0:
                    self.model.change_reaction_bounds(substrate_uptake_id,
                                                  lower_bound = 0, upper_bound = rate)
                else:
                    self.model.change_reaction_bounds(substrate_uptake_id,
                                                          lower_bound=rate, upper_bound=0)
                self.model.optimize()

                #if the model is not optimal revert changes and continue
                if self.model.solver.status != 'optimal':
                    print('infeasible for substrate id', substrate_uptake_id)
                    # make sure infeasibility decreases the error
                    fluxes_df = pd.DataFrame(columns=fluxes_df.columns)
                # calculate fitness (sum of simulation error to reactions with data)
                else:
                    for rxn_id in fluxes_df.columns[1:]:
                        if rxn_id in self.model.reactions:
                            fluxes_df.iloc[-1, fluxes_df.columns.get_loc(rxn_id)] = self.model.reactions.get_by_id(rxn_id).flux

            self._change_kcat_values_for_individual(individual, revert=True)
            error += [self._calculate_simulation_error(fluxes_df, substrate_uptake_id)]

            # reset substrate_uptake_rate
            if self.substrate_uptake_rates[substrate_uptake_id][0]<0:
                self.model.change_reaction_bounds(substrate_uptake_id,
                                              lower_bound=0, upper_bound=1e6)
            else:
                self.model.change_reaction_bounds(substrate_uptake_id,
                                              lower_bound=-1e6, upper_bound=0)
        #average fitness:
        fitness = float(np.nanmean(error))
        individual.r_squared = fitness
        individual.fitness.values = [fitness]


        #revert kcat_changes
        self._change_kcat_values_for_individual(individual, revert = True)

        param_fit = self.init_fitness()
        new_fitness = MyFitness()#weights = param_fit["weights"])
        new_fitness.weights = param_fit["weights"]
        new_fitness.values = tuple([float(fitness)])
        return new_fitness
    
    
    
    ##########################################################################
    # CUSTOM HELPER FUNCTIONS
    ##########################################################################

    def _reset_translational_sector(self):
        self.model.constraints[self.model.TOTAL_PROTEIN_CONSTRAINT_ID].ub = self.totprot_start_ub
        self.model.sectors.get_by_id('TranslationalProteinSector').intercept = self.tps_intercept

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


    def _change_kcat_values_for_individual(self, individual,
                                           kcat_values:list = None,
                                           revert:bool=False):
        if kcat_values is None:
            kcat_values = individual.kcat_list

        enz_to_evaluate = self.enzymes_to_evaluate
        if not revert:
            enz_to_evaluate = {enz_id:
                {
                'reaction':self.enzymes_to_evaluate[enz_id]['reaction'],
                'kcats':{individual.directions[i]:kcat_values[i]}
                } for i, enz_id in enumerate(individual.enzymes_to_eval)
            }
        for enz_id, enz_info in enz_to_evaluate.items():
            self._change_kcat_value_in_model(enz_info['reaction'],
                                             enz_id,
                                             list(enz_info['kcats'].keys())[0],
                                             list(enz_info['kcats'].values())[0])

    def _change_kcat_value_in_model(self, rxn_id:str,
                                    enzyme_id:str,
                                    direction:Literal['f', 'b'],
                                    kcat:float) -> None:

        if enzyme_id not in rxn_id:
            rxn_id = f"CE_{rxn_id}_{enzyme_id}"
        rxn = self.model.reactions.get_by_id(rxn_id)
        if direction == 'b':
            var = rxn.reverse_variable
        else:
            var = rxn.forward_variable
        self.model.constraints[f'EC_{enzyme_id}_{direction}'].set_linear_coefficients({
            var: (1/kcat)})


    def _calculate_simulation_error(self, flux_df: pd.DataFrame,
                                    substrate_reaction:str):
        error = []
        weights = []

        if len(flux_df)==0:#if all rates were infeasible: r_squared should be really bad
            return -1
        for rxn in self.reactions_with_data[substrate_reaction]:
            #only select the rows which are filled with data
            ref_data_rxn = self.valid_data[substrate_reaction].dropna(axis = 0, subset = rxn)

            #if there are no reference data points, continue to the next reaction
            if len(ref_data_rxn) == 0: continue
            r_squared = calculate_r_squared_for_reaction(reaction_id = rxn,
                                                         validation_data = ref_data_rxn,
                                                         substrate_uptake_id = substrate_reaction,
                                                         fluxes = flux_df)
            error += [r_squared]
            if not np.isnan(r_squared):
                # error += [r_squared]
                if rxn in self.weights.keys(): weights.append(self.weights[rxn])
                else: weights.append(1)
        if len(error) == 0: return np.NaN
        np.average(error, weights=weights)
        return np.average(error, weights = weights)