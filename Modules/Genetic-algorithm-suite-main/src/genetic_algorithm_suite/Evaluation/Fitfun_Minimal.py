"""
Minimal example for a fitness function class

Objective: Maximize the number of reaction deletions supporting a minimal growth rate

"""

import random
from pathlib import Path
from os.path import dirname, abspath

# set standard paths
FILE_PATH = Path(abspath(dirname(__file__)))
DATA_PATH = FILE_PATH.parents[0].joinpath("Data")

# seed random number generator
random.seed()

class FitnessEvaluation():
    
    def __init__(self, model=None, fixed_attr_list=[], processes=2):
        """Initialize fitness evaluation class for a genetic algorithm
        
        Inputs:
            :param cobra.core.Model model: Metabolic model in COBRA format
            :param list fixed_attr_list: Identifiers of attributes not to be used as solution variables
            :param int processes: Number of workers available (unused here)
        
        """
                    
        # save list of fixed attributes
        self.fixed_attr_list = fixed_attr_list
        
        # initialize metabolic model
        self.init_model(model)
        
        
    ##########################################################################
    # NECESSARY CUSTOM FUNCTION
    ##########################################################################
    
    def init_model(self, model=None):
        """Pass and save the metabolic model
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
        
        # example property: number of deleted reactions related to a gene
        for individual in pop:
            # determine deleted reactions
            reactions_deleted = self._determine_deleted_reactions(individual)
            # determine deleted reactions associated with a gene
            gene_related_reactions = [
                r
                for r in reactions_deleted
                if len(self.model.reactions.get_by_id(r).genes)>0
                ]
            # save property
            ind_props.append({
                'number_gene_related_reaction_deletions': len(gene_related_reactions)
                })
        
        return ind_props
    
    
    ##########################################################################
    # DEAP-RELATED FUNCTIONS
    ##########################################################################
    def attribute_generator(self, probability=0) -> int:
        """Generates an attribute of a DEAP individual and is part of the 
        mutation operator
        
        Inputs:
            :param float probability: probability for returning a "0" as attr
        
        Outputs:
            :param int attr: representation of an attribute
        
        """
        
        if random.random() >= probability:
            attr = 1
        else:
            attr = 0
        
        return attr
    
    
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
        
    
    def init_attribute(self) -> dict:
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
            
        Outputs:
            :param dict param_attr: 
        
        """
        
        param_attr = {}
        
        # save attribute details for matching attributes to its model identity
        # here each model reaction is an attribute (solution variable)
        self.individual_attr_list = [
            {'id': r.id, 'type': 'deletion'}
            for r in self.model.reactions
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
        

        # create individuals representing metabolic genes as variables/targets
        # creator.create("Individual", list, fitness=creator.FitnessMin)
        param_ind["individual_type"] = list
        
        return param_ind
            
    
    
    ##########################################################################
    # FITNESS FUNCTION
    ##########################################################################
    def eval_fitness(self, individual) -> int:
        """Evaluate the fitness of an individual. It returns a tuple of 
        objective values.
        
        As a simple example, the objective is to maximize the number of 
        deleted reactions not associated with a gene. Any model simulations may
        be performed here to compute the objective function values.
        
        Minimize the computational workload in this function, since it will be
        frequently evaluated
        
        Inputs:
            :param list individual: a list of variables representing an individual (solution)
            
        Outputs:
            :param int fitness: fitness of the current individual evaluated
        
        """

        # Here, the simple objective is to maximize the number of deleted reaction
        # not associated with a gene while allowing for a reasonable growth on glucose

        # apply knockouts and compute metabolic functionalities
        
        with self.model as model:
            # determine deleted reactions
            reactions_deleted = self._determine_deleted_reactions(individual)
            
            # set reaction knockouts in the model
            for r in reactions_deleted:
                model.reactions.get_by_id(r).knock_out()
            
            # simulate the growth rate
            mu = model.slim_optimize(error_value=0)
  

            # calculate fitness (or objective function values)
            # guarantee a minimal growth rate, otherwise set the fitness to zero
            fitness = len(reactions_deleted)*(mu>0.3)
            
        
        # return a tuple of one element
        return tuple([float(fitness)])
    
    
    
    ##########################################################################
    # CUSTOM HELPER FUNCTIONS
    ##########################################################################
    
    def _determine_deleted_reactions(self, individual) -> list:
        """Determine deleted reactions from the set of attributes of an individual
        
        Inputs:
            :param list individual: attributes of an individual (solution)
         
        Outputs:
            :param list reactions_deleted: Identifier of deleted reactions
            
        """
        
        reactions_deleted = [
                attr['id']
                for attr, attr_stat in zip(self.individual_attr_list, individual)
                if attr_stat==0
                ]
        
        return reactions_deleted
        
        
        
