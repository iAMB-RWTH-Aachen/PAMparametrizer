import os
script_directory = os.path.dirname(os.path.abspath(__file__))
cwd = os.getcwd()


os.chdir(os.path.dirname(os.path.dirname(os.path.dirname(script_directory))))
from Modules.genetic_algorithm_parametrization import GAPOGaussian as GAPOGauss
from Modules.genetic_algorithm_parametrization import GAPOUniform
from Scripts.pam_generation import setup_toy_pam

os.chdir(cwd)