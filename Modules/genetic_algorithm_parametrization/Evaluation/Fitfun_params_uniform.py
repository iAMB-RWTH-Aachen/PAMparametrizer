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

from .Fitfun_params_gaussian import FitnessEvaluation

# set standard paths
FILE_PATH = Path(abspath(dirname(__file__)))
DATA_PATH = FILE_PATH.parents[0].joinpath("Data")
DIFUSSIONLIMIT = 1e6

# seed random number generator
random.seed()

class FitnessEvaluation(FitnessEvaluation):



    def __init__(self, *args, **kwargs):

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

        super().__init__(**kwargs)
        
    ##########################################################################
    # NECESSARY CUSTOM FUNCTION
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
            kcats = [kcat
                     if kcat<DIFUSSIONLIMIT else DIFUSSIONLIMIT
                     for kcat in np.random.lognormal(mean = np.log10(self.KCAT_MU),
                                                     sigma = np.log10(self.KCAT_SIGMA),
                                                     size = self.NUM_KCATS)
                     ]
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

    ##########################################################################
    # CUSTOM HELPER FUNCTIONS
    ##########################################################################
        
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
        if toolbox is not None:
            # mutate an individual with a mutation rate based on the sensitivity of the individual enzymes
            # the new value is samples from a gaussian distribution with mu being the original kcat value and
            # sigma being related to the kcat value to stay in sync with the order of magnitude of the original kcat
            return toolbox.mutate([kcat])[0]
        else:
            # scale needs to be positive, so prevent negative values
            # loc = mean, scale = sd, sd is defined as kcat/10 to make sure sd is in the same order of magntitude
            # as the kcat is
            max_kcat = kcat*2 if kcat*2<DIFUSSIONLIMIT else DIFUSSIONLIMIT
            return float(np.random.uniform(0,max_kcat))

    def _mut_kcat_uniform(self, kcat_list:list, indpb:float):
        """
        Mutate kcat list by sampling from a uniform distribution ranging from 0 to 2*kcat
        Mutation probability is determined by indpb

        :param kcat_list: list of kcats to mutate
        :param indpb: Independent probability for each attribute to be mutated.
        :return: new_kcats: list with mutated kcat values
        """
        if random.random() < indpb:
            new_kcats = []
            for kcat in kcat_list:
                max_kcat = kcat*2 if kcat * 2 < DIFUSSIONLIMIT else DIFUSSIONLIMIT
                new_kcats.append(np.random.uniform(0, max_kcat))
        else:
            new_kcats = kcat_list
        return new_kcats