# Workflow for Protein Allocation Model Parameter Estimation
This subdirectory contains the framework which performs the parametrization workflow. The workflow starts
from an initial parametrization. Based on model simulations and the effect of individual enzymes on the 
simulated growth rate, enzymes are selected. The parameters of these enzymes will be optimized by a genetic
algorithm, which aims to maximize the fit of the simulations to experimentally observed phenotypes by means
of an R<sup>2</sup> value. The whole workflow is repeated until a user-defined maximum number of iterations or error 
threshold is reached.

## Software Structure
The PAM parametrizer is organized as follows:
1. `PAM_parametrizer`: contains the `PAMParametrizer` object with all the functionalities
2. `PAM_data_classes`: contains all data classes to store metadata, hyperparameters and results
   * `ValidationData`: experimental data used to determine the R^2 value and the reactions which should be used for validation/plotting
   * `Hyperparameter`: all hyperparameters for the parametrizer and the genetic algorithm
   * `ParametrizationResults`: results of the parametrization workflows, such as fluxes, sensitivities and best parameter sets
   * `FluxResults`: class storing the simulation results (fluxes and error to experimental measurements) for a single substrate uptake reaction. Stored in ParameterizationResults.flux_results.

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
>Save the resulting png file and parameters to excel


## Hyperparameters and other inputs
