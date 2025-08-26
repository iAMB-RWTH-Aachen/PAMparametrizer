"""
Minimal example for a fitness function class

Objective: Maximize the number of reaction deletions supporting a minimal growth rate

"""

import random

import deap.base
import numpy as np
from typing import List, Optional

from .Fitfun_params_gaussian import FitnessEvaluation, DIFUSSIONLIMIT

# seed random number generator
random.seed()

class FitnessEvaluation(FitnessEvaluation):



    def __init__(self, *args, **kwargs):
        super().__init__(**kwargs)
        
    ##########################################################################
    # NECESSARY CUSTOM FUNCTION
    ##########################################################################

    def attribute_generator(self, probability: float=0,
                            kcat_list:Optional[List[float]] = [],
                            kcat_bounds: Optional[List[dict]] = {}
                            ) -> int:
        """Generates an attribute of a DEAP individual and is part of the 
        mutation operator
        
        Args:
            robability (float): probability for returning a "0" as attr
            kcat_list (List[float]): list of kcat values
        
        Returns:
            attr (int): representation of an attribute
        
        """

        if kcat_list == []:
            kcats = [kcat
                     if 1/kcat<DIFUSSIONLIMIT else 1/DIFUSSIONLIMIT
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
                kcat_constraint = kcat_bounds[i]
                new_kcat = self._mutate_kcat_value(kcat,**kcat_constraint)
                kcat_list[i] = abs(new_kcat) #new kcat value should always be positive

        return kcat_list

    ##########################################################################
    # CUSTOM HELPER FUNCTIONS
    ##########################################################################
        
    def _mutate_kcat_value(self, kcat:float,
                           min_kcat: float = 1e-6,
                           max_kcat:float = DIFUSSIONLIMIT,
                           sensitivity:float = 0.5,
                           toolbox:deap.base.Toolbox = None) -> float:
        """
        Helper function to mutate a kcat value based on a user-defined uniform scale.
        Can be used both with the deap toolbox as well as without. In the latter case,
        a new kcat value will be sampled from np.random.uniform distribution.

        Args:
            kcat (float): kcat value to mutate
            min_kcat (Optional float): minimal kcat value, defaults to 1e-6 (not 0 to prevent arithmetic errors)
            max_kcat (Optional float): maximal kcat value, defaults to diffusion limit (1e6 1/s)
            sensitivity (Optional float): sensitivity of the enzyme related to the kcat to determine mutation probability
            toolbox (Optional deap.base.toolbox): deap toolbox

        Returns
            new_kcat (float): mutated kcat value (sampled from uniform distribution)
        """
        if toolbox is not None:
            # mutate an individual with a mutation rate based on the sensitivity of the individual enzymes
            # the new value is samples from a gaussian distribution with mu being the original kcat value and
            # sigma being related to the kcat value to stay in sync with the order of magnitude of the original kcat
            return toolbox.mutate([kcat], min_kcat = min_kcat, max_kcat = max_kcat)[0]
        else:
            # scale needs to be positive, so prevent negative values
            # loc = mean, scale = sd, sd is defined as kcat/10 to make sure sd is in the same order of magntitude
            # as the kcat is
            max_kcat = kcat * 2 if 1 / (kcat * 2) < max_kcat else 1 / max_kcat
            return float(np.random.uniform(max_kcat, 1/min_kcat))

    def _mut_kcat_uniform(self,
                          kcat_list:List[float],
                          indpb:float,
                          min_kcat: float = 1e-6,
                          max_kcat:float = DIFUSSIONLIMIT) -> List[float]:
        """
        Mutate kcat list by sampling from a uniform distribution ranging from 0 to 2*kcat
        Mutation probability is determined by indpb

        Args:
            kcat_list List[float]: list of kcats to mutate
            indpb (float): Independent probability for each attribute to be mutated.
            min_kcat (Optional float): minimal kcat value, defaults to 1e-6 (not 0 to prevent arithmetic errors)
            max_kcat (Optional float): maximal kcat value, defaults to diffusion limit (1e6 1/s)

        Returns
            new_kcats (list[float]): list with mutated kcat values
        """
        if random.random() < indpb:
            new_kcats = []
            for kcat in kcat_list:
                max_kcat = kcat*2 if 1/(kcat * 2) < max_kcat else 1/max_kcat
                new_kcats.append(np.random.uniform(max_kcat, 1/min_kcat))
        else:
            new_kcats = kcat_list
        return new_kcats