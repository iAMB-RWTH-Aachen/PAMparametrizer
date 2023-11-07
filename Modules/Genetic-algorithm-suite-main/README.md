# Genetic Algorithm for Metabolic Model Optimization
Framework for using a genetic algorithm to optimize problems involving simulations of metabolic models in COBRA format. For more information about the terminology and principles of genetic algorithms for metabolic optimization refer to this [paper](https://doi.org/10.3390/metabo8020033).  
The [DEAP toolbox](https://github.com/DEAP/deap) ([DEAP documentation](https://deap.readthedocs.io/en/master/index.html)) is used for basic genetic algorithm functionalities. The initial usecase of this framework was to compute a minimal set of metabolic genes preserving certain core metabolic functionalities.

## Installation 
The Genetic Algorithm Suite requires an installation of the Gurobi solver and a valid (academic) licence. Refer to the [Gurobi website](https://www.gurobi.com/) for further instructions.

It is recommended to create a new virtual environment (conda -> tested, virtualenv) with Python 3.8 or greater for installing the Genetic Algorithm Suite. For creating a conda environment run `conda create --name env_ga python=3.9` from a shell, and activate it with `conda activate env_ga`.
1. Clone, fork, or download this repository
2. Browse to the main directory of this repository
3. Run `python setup.py develop` or `python setup.py install`

cobrapy and gurobipy for handling metabolic models, as well as deap are automatically installed. If problems occur when running the Genetic Algorithm Suite, try to separately install specififc versions of these packages in a new environment before installing the Genetic Algorithm Suite.
1. Run `pip install cobrapy==0.25.0`
2. Run `pip install gurobipy==9.5.2`. The gurobipy package needs to support your installed version of the gurobi solver.
3. Run `pip install deap==1.3.1`

## Use
Run the example "minimal_example.py" in "/examples" (`python minimal_example.py` from a shell) to test the Genetic Algorithm Suite and get an overview of its functionalities. The corresponding fitness function class "Fitfun_Minimal" can be found in "src/genetic_algorithm_suite/Evaluation". This is a minimal working example and shows the structure and all requirements demanded from a custom fitness function class by the Suite's core routines. "Fitfun_Minimal.py" may be copied, edited, and extended for solving any metabolic model-related problem beyond metabolic reduction.

The GAMO class (Genetic Algorithm for Metabolic Optimization) is the core of the Genetic Algorithm Suite and initializes an optimization problem. Refer to the example scripts for information about the various keyword arguments needed to be passed for setting up the genetic algorithm.

The output of the Genetic Algorithm Suite is stored in the directory specified by the keyword argument "folderpath_save" in the GAMO class and consists of three files:
1. .pickle with the initialized core module
2. .json with metadata about the optimization run and the latest population. This file is used to restart the optimization, eventually with altered optimization parameters, and the path to this file is passed to the "restart" method.
3. An Excel file including human-readable information about the optimization run and the final population

If intermediate results are not overwritten ("overwrite_intermediate_results=False"), this output is generated for each gene-flow-event, i.e. after each sequence of multiple generations distributed to multiple parallel workers.

For applying the comprehensive minimal metabolic genome fitness function class refer to "metabolic_genes_reduction.py" in "/examples". The class itself can be found in "src\genetic_algorithm_suite\Evaluation\Metabolic_Genes_Reduction.py" and again may be copied and edited for using a metabolic model other than [iML1515](http://bigg.ucsd.edu/models/iML1515) or changing the enforced metabolic functionalities.



