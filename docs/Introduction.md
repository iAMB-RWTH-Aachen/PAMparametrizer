# PAMparametrizer - parametrizing Protein Allocation Models with Flux Data
*******

## Why do we need a PAMparametrizer?
The Protein Allocation Model is a constraint-based model which includes constraints connecting enzyme concentrations
to reaction rates, and limits the total enzyme abundance in a system. As a result, the metabolic phenotype can be modeled
with more accurate without considerable increase in computational costs 
([Alter et al. 2021](https://journals.asm.org/doi/10.1128/mSystems.00625-20)). The relation between reaction rates and 
enzyme concentrations depends on the turnover number of the associated enzymes (the so-called k<sub>cat</sub>). Even though there 
are multiple databases available which provide these turnover numbers, e.g. [BRENDA](https://www.brenda-enzymes.org/) 
and [SABIO-RK](https://sabiork.h-its.org/). Nevertheless, with the increasing number of identified enzymes and 
newly created models of organisms, these number of enzymes without an experimentally determined k<sub>cat</sub> value 
increases. Furthermore, the experimentally determined *in vitro* k<sub>cat</sub> values often do not agree with the 
*in vivo* k<sub>cat</sub> values in biological systems. Multiple methods are available to parametrize protein-constrained
models with omics data (Wilken et al. 2022, Davidi et al. 2016) or machine learning methods (Wendering2023, Heckmann2018, Li2022).
The latter methods unfortunately do not perform well an additional fitting step, as they are trained on the available in vitro data.

If you are interested in designing experiments or understanding novel chassis strains with a new PAM, large amounts of quantitative omics 
data are often not available and/or difficult to measure. In these cases, we would like to design experiments with only
limited amount of experimental effort. This is where the PAMparametrizer comes in. The PAMparametrizer is designed to fit
a set of parameters, obtained by machine learning, to simple phenotypic data such as exchange rates and growth rates. As
a result PAMs can be used early in the development of novel processes with novel microbes.


## What does the PAMparametrizer do?
The PAMparametrizer tries to solve this issue by using just simple exchange rates as input for the parametrization. It makes
use of a powerful genetic algorithm to mutate the initial parameter set until the fit of the model simulations to experimental
measurements is maximized. In other words: it performs 'evolution on the computer'. It takes the following steps:

1. Simulating phenotypes
2. Select the enzymes with the largest effect on the growth objective using sEnz
3. Optimize the k<sub>cat</sub> values of the selected enzymes using a genetic algorithm
4. Evaluate and save the final results
5. Repeat 1-5 until a maximum fitness or maximum number of iterations is reached

## What can you find where in this repository?
This repository contains not only the source code, but also examples and scripts which were used in **INSERT PUBLICATION HERE**.

### Data and Scripts for van den Bogaard, et al (2025)

### PAMparametrizer functionality
Besides the code which are used to perform the analyses as described in ....., this repository also contains the 
source code and other useful scripts for anyone who wants to create a PAM.

- **Modules** 
    - *genetic_algorithm_parametrization*: package to run the [genetic algorithm](/Genetic_algorithm).
    - *PAM_parametrizer*: package to run the [parametrization workflow](/PAM_param)
    - *util*: some utility functions which occur both in the genetic algorithm as the PAMparametrizer, such as determining calculating and R<sup>2</sup> and parametrizing the translational sector
- **Models**: models to use as test cases
- **Scripts**
    - *Protein_sectors*: scripts with examples on how to parametrize the protein sectors
    - *Testing*: example scripts on how to set up and run the parametrization workflow for different models
    - *pam_generation_uniprot_id.py*: functions to easily set up a PAM
    - *parse_kcat_values_GotEnzymes.ipynb*: jupyter notebook explaining how to do the initial parametrization with the GotEnzymes database
- **tests**: all unit tests which the framework should pass. This can also be used to see how the algorithm should be used and should behave.

IMPORTANT: all code should be run from the main repository via the terminal. E.g. `python -m Scripts.Testing.pam_parametrizer_iML1515.py`

## Code structure
For a more detailed look on the interdependencies of the framework and the structure of the packages, please have a look at
the UML diagram of the entire software structure:

![PAMparamUML](PAMparametrizer_UML.svg)

## Dependencies
Both the PAMparametrizer package as the genetic algorithm use the modeling framework [PAModelpy](https://github.com/iAMB-RWTH-Aachen/PAModelpy)
to handle the Protein Allocation Models. This package is able to incorporate complex gene-protein-reaction associations, 
such as AND and OR relations, incorporating the complexity present in biological systems. This package also includes 
the functionalities to perform the enzyme sensitivity analysis sEnz (refer to senz paper).

The genetic algorithm makes use of the [DEAP toolbox](https://github.com/DEAP/deap). This modular toolbox includes all
the tools needed to create different types of evolutionary algorithms, such as the genetic algorithm.

The full list of dependencies is as follows:
- PAModelpy==0.0.4.8
- matplotlib==3.8.4
- deap==1.4.1y
- random
- timepip
- gurobipy==11.0.1
- scipy==1.13.0
- pytest==8.2.0
- seaborn==0.13.2
- scikit-learn==1.4.2
- plotly

All dependencies can be installed in one go by downloading this repository and running:

`pip install -r requirements.txt`

from the main directory

## License
Copyright institute of Applied Microbiology, RWTH Aachen University, Aachen, Germany (2023)

PAModelpy is free of charge open source software, which can be used and modified for your particular purpose under the [MIT](https://opensource.org/license/mit/)
or [Apache 2.0](https://www.apache.org/licenses/LICENSE-2.0) of the users choice.

Please note that according to these licenses, the software is provided 'as is', WITHOUT WARRANTY OF ANY KIND, 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.