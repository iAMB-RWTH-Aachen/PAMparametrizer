"""
Genetic algorithm (GA) for the prediction genome reduction paths based on metabolic models
- DEAP (https://github.com/DEAP/deap) is used as a base GA
- GA can be generalized to any optimization problem incorporating M-models

"""
import numpy as np
# Disable gurobi logging output
from .core_parametrization_gaussian import GAPO
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


class GAPOUniform(GAPO):
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _init_deap_toolbox_mutation(self, toolbox):
        # only change the mutation function to uniform sampling
        toolbox.register("mutate", self.FitEval._mut_kcat_uniform ,indpb=self.mutation_rate)

        return toolbox
