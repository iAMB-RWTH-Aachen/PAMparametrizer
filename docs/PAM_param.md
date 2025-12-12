# Workflow for Protein Allocation Model Parameter Estimation
This subdirectory contains the framework which performs the parametrization workflow. The workflow starts
from an initial parametrization. Based on model simulations and the effect of individual enzymes on the 
simulated growth rate, enzymes are selected. The parameters of these enzymes will be optimized by a genetic
algorithm, which aims to maximize the fit of the simulations to experimentally observed phenotypes by means
of an R<sup>2</sup> value. The whole workflow is repeated until a user-defined maximum number of iterations or error 
threshold is reached. As a final step, the amount of unused enzyme sector at zero growth is optimized to determine the 
quantitative protein burden.

## Software Structure
The PAM parametrizer is organized as follows:
1. `PAM_parametrizer`: contains the `PAMParametrizer` object with all the functionalities
2. `PAM_data_classes`: contains all data classes to store metadata, hyperparameters and results
   * `ValidationData`: experimental data used to determine the R<sup>2</sup> value and the reactions which should be used for validation/plotting
   * `Hyperparameter`: all hyperparameters for the parametrizer and the genetic algorithm
   * `ParametrizationResults`: results of the parametrization workflows, such as fluxes, sensitivities and best parameter sets
   * `FluxResults`: class storing the simulation results (fluxes and error to experimental measurements) for a single substrate uptake reaction. Stored in ParameterizationResults.flux_results.
3. `KcatConstraintConfig`

The parametrization workflow can make use of experimental measurements of growth on different substrates.
This is enabled by providing a separate ValidationData instance for each substrate. This way, each substrate 
is stored together with the corresponding experimental data and reactions which should be used for validation.
The user can define one *main* substrate which is purely used to visualize the progress of the parametrization in an
effective manner. For each substrate, a separate FluxResults object will be created and stored in ParametrizationResults.
As a consequence, the simulations and corresponding errors can be stored for each substrate individually.

## PAM parametrizer workflow
The `PAMparametrizer` is run with the `run` function. The workflow is organized as follows (in pseudocode):

> For each carbon source:
>> Calculate relation between substrate uptake rate and translational protein sector
> 
> Plot validation data  
> 
> if number_experimental_measurement > hyperparameter.bin_resolution: Sample the validation data  
> 
> if binned == before: Perform iteration in bins  
> 
> while (iteration <= hyperparameters.threshold_iteration) & (final_error <= hyperparameters.threshold_error):
>>0. Initialize result objects in `ParametrizationResults`.
>>1. Run simulations over range of substrate uptake rates
>>2. Based on the sensitivity analysis described in van den Bogaard, et al. (2024), select enzymes to evaluat
>>3. Initiate genetic algorithm with sampled validation data and associated substrate uptake rates
>>4. genetic_algorithm.start() (more details in Genetic Algorithm README file)
>>5. Parse the results of the genetic algorithm, save intermediate results to excel file
>>6. Reparametrize the model with optimized parameters
>>7. Plot the progress
>Optimize the unused enzyme sector at zero growth
>Save the resulting png file and parameters to excel


## Data objects
The PAMparametrizer is a complex framework requiring different types of input from the user and generates different types
of output. To orchestrate the data movement through the framework, the inputs and outputs of the PAMparametrizer
are organized in four different data classes in the `PAM_data_classes.py` file. Below follows a description of each of 
the classes, and how it is used by the framework.

### SectorConfig <USER INPUT>
The main power of the PAM is that it contains sectors which represent proteins which are not included explicitely in 
protein-constrained metabolic models, but do occupy a variable amount of protein space. Normally, 2 sectors are added:
the UnusedEnzymeSector and the TranslationalProteinSector. The amount of proteins in these sectors are related to the (carbon)
substrate uptake rate. This poses a challenge: the relation between carbon substrate and the amount of proteins in these 
sectors is unique for each substrate. Nevertheless, when a metabolic model is not limited by the amount of proteins,
there is a linear relation between the carbon uptake and the growth rate, depending on the way the microbe metabolizes the substrate.
We can thus use the relationship between the growth rate and protein sector in not protein-limited regimes to compute the 
relation between the protein sector and the various substrates. This relationship is stored in the input Excel file for the
PAMparametrizer and provided to the PAMparametrizer object using the SectorConfig TypedDict:

- **sectorname**: The sector identifier as present in the PAM
- **slope**: relationship between sector (g/gCDW) and the growth rate (1/h)
- **intercept**: sector (g/gCDW) at zero growth
- **substrate_range**: range of substrate uptake rates in which the relation between sector and growth rate is valid

### ValidationData <USER INPUT>
The PAMparametrizer needs exchange fluxes to determine how good the simulations with the new parameter set are. These flux rates
can be measured in different conditions, most commonly using different carbon sources. This poses a complex problem:
how do we provide the framework with not only the experimental measurements, but also the precise conditions in which the
measurements have been performed, and which reactions are most representative to see the progress? The ValidationData object 
tries to streamline this data handling issue. For each condition to which the parameter set should be optimized, the corresponding 
ValidationData object should be created. It stores the following information:

- **id**: substrate uptake identifier
- **valid_data**: dataframe containing the experimental measurements
- **validation_range**: range of substrate uptake range to consider for validation
- **reactions_to_validate**: reactions to use for error calculation
- **reactions_to_plot**: reactions to use for visualization
- **sector_config**: configuration of all sectors related to the substrate uptake rate (as SectorConfig)

When the sector_configs is left empty, it will be automatically determined using the relation between the sector and 
the growth rate. The PAMparametrizer recalculates the sector parameters using the relation between the growth rate and 
substrate uptake rate of a PAM with a large protein pool (similar to a genome-scale model). If the relation between the 
sector parameters and the growth rate are unknown, only the translational sector is assumed to be related to the substrate uptake rate.
The sector parameters are then obtained, assuming a similar growth rate to translational protein content relationship 
as established for *E. coli*.

In the PAMparametrizer object, all the ValidationData objects are stored and used by the framework to validate the new
parametersets.

### HyperParameters <USER INPUT>
A genetic algorithm is a very flexible tool for optimizing a set of parameters based on some type of error calculation.
It can, however, be sensitive to different hyperparameters, such as the mutation probability and rate, the number of
generations, the population size, etc. On top of that, the PAMparametrizer itself also needs some manual input. For
examples the number of iterations, the convergence threshold, the number of kcats to optimize, etc. For easy handling 
and easy reproduction, all this information is saved in the HyperParameters data object. 

Below follows a summary off all the inputs you can change in this object, including an explanation of its effects.


| Attribute                      | Default value            | Explanation                                                                                                                                                                                   |
|--------------------------------|--------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **PAMparametrizer**            |
| threshold_error                | 0.9                      | R<sup>2</sup> value at which workflow stops                                                                                                                                                   |
| threshold_iteration            | 100                      | max number of iterations                                                                                                                                                                      |
| threshold_convergence          | 1e-3                     | max difference between subsequent errors to consider the error to be converging                                                                                                               |
| threshold_nmbr_convergence     | 3                        | number of final error to consider for convergence                                                                                                                                             |
| number_of_kcats_to_mutate      | 5                        | number of k<sub></sub> to optimize per iteration                                                                                                                                              |
| genetic_algorithm_file_base    | 'genetic_algorithm_run_' | base of files in which the genetic algorithm results are saves                                                                                                                                |
| filename_extension             | ''                       | filename extension to recognize output files with                                                                                                                                             |
| genetic_algorithm              | GAPOUniform              | personalized genetic algorithm object <br/>(if you want to change something in the genetic algorithm structure)                                                                               |
| genetic_algorithm_hyperparams  | dict()<br/> see below    | dictionary containing all hyperparameters for the genetic algorithm<br/>see below for all inputs                                                                                              |
| **Genetic Algorithm**          |
| mutation_probability           | 0.5                      | probability with which an individual (solution) is mutated in a generation                                                                                                                    |
| mutation_rate                  | 0.5                      | probability with which an attribute (e.g. gene) of an individual is mutated                                                                                                                   |
| population_size                | 10                       | number of individuals (solution) per population                                                                                                                                               |
| crossover_probability          | 0.8                      | probability with which two indivduals/offsprings are crossed over                                                                                                                             |
| number_generations             | 20                       | number of consecutive generations per gene flow event                                                                                                                                         |
| number_gene_flow_events        | 2                        | number of gene flow events, <br/>i.e. merging of multiple populations independently evolved on parallel workers                                                                               |
| init_attribute_probability     | 0                        | probability with which the initial population is mutated                                                                                                                                      |
| fitness_class                  | 'Fitfun_params_uniform'  | name of the file containing the personalized fitness function                                                                                                                                 |
| processes                      | 2                        | number of parallel workers                                                                                                                                                                    |
| time_limit                     | 600                      | time limit in seconds                                                                                                                                                                         |
| error_weights                  | dict()                   | reaction which should have a different impact on the error calculation <br/>than other reactions<br/>e.g. {'EX_ac_e':5}:acetate has 5 times more impact on the error than the other reactions |
| folderpath_save                | Path(r"Results")         | path for saving results                                                                                                                                                                       |
| overwrite_intermediate_results | True                     | if true, saved intermediate results are overwritten                                                                                                                                                                                              |
| print_progress                 | True                     | if True, progress of the genetic algorithm is printed                                                                                                                                                                                              |


### ParametrizationResults <Intermediate OUTPUT>
The PAMparametrizer does not only need a complex input, but also generated a complex output. For this reason, the output
is organized in two different data structures: (1) ParametrizationResults, to save the (intermediate) results of the entire
framework, and (2) FluxResults, to save the (intermediate) results of the simulation of various conditions.

The ParametrizationResults object is initiated upon initiation of the PAMparametrizer. It saves the following information:
- **substrate_uptake_reactions**: list of all substrate uptake reactions to simulate
- **esc_df**: dataframe to store the enzyme sensitivities
- **sensitive_enzymes**: dataframe to store the most sensitive enzymes selected for each iteration
- **flux_results**: dictlist to store the FluxResults objects for all conditions
- **best_individuals**: dataframe to store the best individuals for genetic algorithm optimization
- **computational_time**: dataframe to store the computational costs of each run
- **final_errors**: dataframe to store the final error of each simulation

### FluxResults <Intermediate OUTPUT>
The FluxResults object is stored in the ParametrizationResults. It saves the results, mostly flux rates, of the simulations
of a specific conditions, which can be used at a later time point to calculate the error to experimental data (stored in
ValidationData). It stores the following:

- **id**: substrate uptake identifier
- **error_df**: dataframe to store the error for this specific conditions
- **fluxes_df**: dataframe to store the reaction fluxes required to calculate the error (defined by `ValidationData.reactions_to_validate`)
- **substrate_range**: a list with the substrate uptake rates to validate. The substrate uptake rates are obtained from `ValidationData.sampled_valid_data`