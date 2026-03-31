# PAMparametrizer - parametrizing Protein Allocation Models with Flux Data
*******

## Why do we need a PAMparametrizer?
The Protein Allocation Model is a constraint-based model that includes constraints connecting enzyme concentrations
to reaction rates, and limits the total enzyme abundance in a system. As a result, the metabolic phenotype can be modeled
with more accurate without considerable increase in computational costs 
([Alter et al. 2021](https://journals.asm.org/doi/10.1128/mSystems.00625-20)). The relation between reaction rates and 
enzyme concentrations depends on the turnover number of the associated enzymes (the so-called k<sub>cat</sub>). There 
are multiple databases available which provide these turnover numbers obtained from experimental *in vitro* measurements,
e.g. [BRENDA](https://www.brenda-enzymes.org/) and [SABIO-RK](https://sabiork.h-its.org/). However, with the increasing number of identified enzymes and 
newly created models of organisms, the number of enzymes without an experimentally determined k<sub>cat</sub> value 
increases. Furthermore, the experimentally determined *in vitro* k<sub>cat</sub> values often do not agree with the 
*in vivo* k<sub>cat</sub> values in biological systems. Multiple methods are available to parametrize protein-constrained
models with omics data ([Wilken et al. 2022](https://www.sciencedirect.com/science/article/pii/S1096717622001173), 
[Davidi et al. 2016](https://www.pnas.org/doi/10.1073/pnas.1514240113)) or machine learning methods (
[Wendering et al. 2023](https://www.nature.com/articles/s41467-023-37151-2), 
[Heckmann et al. 2018](https://www.doi.org/10.1073/pnas.2001562117), 
[Li et al. 2022](https://www.nature.com/articles/s41929-022-00798-z)). The former methods require large amounts of 
proteomics and flux measurements. Also, the latter methods unfortunately do not perform well for the parametrization of metabolic models. 
As they are trained on the available *in vitro* k<sub>cat</sub> values, they do need an additional fitting step.

If you are interested in designing experiments or understanding novel chassis strains with a new PAM, large amounts of quantitative omics 
data are often not available and/or difficult to measure. In these cases, we would like to design experiments with a
limited amount of experimental effort. This is where the PAMparametrizer comes in. The PAMparametrizer is designed to fit
a set of parameters, obtained by machine learning, to simple phenotypic data such as exchange rates and growth rates. As
a result, PAMs can be used early in the development of novel processes with novel microbes.


## What does the PAMparametrizer do?
The PAMparametrizer tries to solve this issue by using just simple exchange rates as input for the parametrization. It makes
use of a powerful genetic algorithm to mutate the initial parameter set until the fit of the model simulations to experimental
measurements is maximized. In other words: it performs 'evolution on the computer'. It takes the following steps:

1. Simulating phenotypes
2. Select the enzymes with the largest effect on the growth objective using [sEnz](https://doi.org/10.1093/bioinformatics/btae691)
3. Optimize the k<sub>cat</sub> values of the selected enzymes using a genetic algorithm
4. Evaluate and save the final results
5. Repeat 1-5 until a maximum fitness or maximum number of iterations is reached
6. Optimize the amount of unused enzyme sector at zero growth

## What can you find where in this repository?
This repository contains not only the source code, but also examples and scripts which were used in **INSERT PUBLICATION HERE**.

### PAMparametrizer functionality
Besides the code which are used to perform the analyses as described in **INSERT PUBLICATION HERE**, this repository also contains the 
source code and other useful scripts for anyone who wants to create a PAM.

- **Modules** 
    - *genetic_algorithm_parametrization*: package to run the [genetic algorithm](/Genetic_algorithm).
    - *PAM_parametrizer*: package to run the [parametrization workflow](/PAM_param)
    - *utils*: some utility functions which occur both in the genetic algorithm as the PAMparametrizer, such as determining calculating and R<sup>2</sup> and parametrizing the translational sector. It also contains utilities to set up the parametrizer and to analyse the results
- **Models**: models to use as test cases
- **Figures**: scripts and code for the figures as presented in **INSERT PUBLICATION HERE**
- **Results**:
    - *1_preprocessing*: initial parametrization files
    - *2_parametrization*: output of the PAMparametrizer
      - diagnostics: Excel files with the changed k<sub>cat</sub>, computation time, and sector parameters
      - progress: prognosis plot showing the simulation results after each iteration of the PAMparametrizer
      - proteinAllocationModel *multi* files: the parameter file as used by the latest parametrization run (e.g., with scaled kcat values)
    - *3_analysis*: analysis of the newly parametrized models
- **Scripts**
    - *i1_preprocessing*: scripts with examples on how to parametrize the protein sectors and prepare the files for the PAMparametrizer
    - *i2_parametrization*: example scripts on how to set up and run the parametrization workflow for different models
    - *i3_analysis*: scripts which help with analysing the results of (multiple) parametrization efforts
    - *Shell*: example shell scripts to run the PAMparametrizer on a high-performance computing cluster with SLURM.
    - *Data_requirements*: scripts to analyse the data-dependency of the PAMparametrizer for parametrizing the iML1515 PAM
- **tests**: testsuite for most of the code in this repository. This can also be used to see how the algorithm should be used and should behave.

IMPORTANT: all code should be run from the main repository via the terminal. E.g. `python -m Scripts.Testing.pam_parametrizer_iML1515.py`

## Code structure
For a more detailed look on the interdependencies of the framework and the structure of the packages, please have a look at
the UML diagram of the entire software structure:

![PAMparamUML](PAMparametrizer_UML.svg)

## Dependencies
Both the PAMparametrizer package as the genetic algorithm use the modeling framework [PAModelpy](https://github.com/iAMB-RWTH-Aachen/PAModelpy)
to handle the Protein Allocation Models. This package is able to incorporate complex gene-protein-reaction associations, 
such as AND and OR relations, incorporating the complexity present in biological systems. This package also includes 
the functionalities to perform the enzyme sensitivity analysis [sEnz](https://doi.org/10.1093/bioinformatics/btae691).

The genetic algorithm makes use of the [DEAP toolbox](https://github.com/DEAP/deap). This modular toolbox includes all
the tools needed to create different types of evolutionary algorithms, such as the genetic algorithm.

The full list of dependencies can be found in requirements.txt
All dependencies can be installed in one go by downloading this repository and running:

`pip install -r requirements.txt`

from the main directory

## License
Copyright institute of Applied Microbiology, RWTH Aachen University, Aachen, Germany (2025)

PAMparametrizer is released under both the GPL and LGPL licenses version 2 or later. 
You may choose which license you choose to use the software under.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
or the GNU Lesser General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
