# PAMparametrizer - parametrizing Protein Allocation Models with Flux Data
*******
The PAMparametrizer tries to simplify the creation of Protein Allocation Models (PAMs) (more info on PAMs see [Alter et al. 2021](https://journals.asm.org/doi/10.1128/mSystems.00625-20),
and [van den Bogaard et al. 2024](https://doi.org/10.1093/bioinformatics/btae691)/[PAModelpy](https://github.com/iAMB-RWTH-Aachen/PAModelpy))
for novel organisms. This helps us to improve our understanding about the metabolic strategies of an organisms with limited
amount of quantitative data available. By combining constraint-based modeling with principles of evolution, the PAMparametrizer
can optimize the parametrization of PAMs using simple phenotypic measurements (nutrient exchanges) as a basis.
It uses a powerful genetic algorithm to mutate the initial parameter set until the fit of the model simulations to experimental
measurements is maximized. 

## Installation instructions
After cloning the repository from git, you can easily install the PAMparametrizer as a python package as follows:

```commandline
pip install ./Modules
```

You can use the `-e` flag during the installation in case you want to modify the source code

```commandline
pip install -e ./Modules
```

## Documentation
Documentation including examples and suggestions on how to perform the parametrization can be found in the [documentation](https://pamparametrizer.readthedocs.io/en/latest/)

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

IMPORTANT: all code should be run from the main repository via the terminal. E.g. `python -m Scripts.i2_parametrization.pam_parametrizer_iML1515.py`

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

