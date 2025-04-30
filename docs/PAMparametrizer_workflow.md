# Parametrizing a PAM: from GEM to PAM
*******

Building a Protein Allocation Model (PAM) from a genome-scale metabolic model (GEM) requires several inputs and design choices
from the user. In this workflow, all the steps, requirements, and inputs are listed. All steps are, where possible,
connected to scripts and notebooks which support the process of going from GEM to PAM. Let's get started!

## 1. The Genome-Scale Metabolic Model Sanity Check
The basis of each PAM is a genome-scale metabolic model of high-quality. To create a PAM out of a GEM, the following requirements
should be met:
- 100% mass balanced
- 100% charge balanced

The following is not required, but makes your life easier:
- Reactions annotated with KEGG ids
- Metabolites annotated with KEGG ids
- Genes annotated with KEGG ids
- \>60% of all reactions should be annotated with gene-reaction associations

Also, make sure you have the following information
- Regular expression to recognize gene-ids (e.g., `r'\b([b|s]\d{4})\b'` for the b1234 locus tags in *E. coli*)
- Potential exceptions in identifier/mapping format (an example is `s0001` which is used for unknown gene)

## 2. Building the protein sectors
You can find scripts to help you build the protein sectors in `Scripts/i1_preprocessing`. The scripts can be run in random
order, although we suggest starting with the ActiveEnzymesSector.

### 2.1 Initiating the ActiveEnzymesSector
**Notebook**: `Scripts/i1_preprocessing/0_parse_kcat_values_GotEnzymes.ipynb` (personlization needed)

**Input**:
 - Genome-scale model
 - Default protein sector parameter file `Data/proteinAllocationModel_EnzymaticData_empty.xlsx`
 - **UniProt** mapping between genes and proteins (instructions in notebook)
 - All k<sub>cat</sub>s for the microorganism from **GotEnzymes** (instructions in notebook)


**Defaults**:
 - Reactions or genes which cannot be mapped to any protein identifier will be given a default id, length and molar mass : `Enzyme_<gene>` of `Enzyme_<reaction>`
 - Proteins which cannot be mapped to a k<sub>cat</sub> value are assigned a default value of 13.7 (BRENDA mean value according to [Bar-Even et al. (2011)](https://doi.org/10.1021/bi2002289))
 - Each protein represents a functional enzyme: enzyme complexes are represented by ids which combine all subunits (`peptide1_peptide2`) and the molar masses and lengths are summed
 - Each row represents an individual enzyme-reaction relation 


**Output**:
 - Excel file saved in `Results/1_preprocessing/proteinAllocationModel_EnzymaticData_<your-model>_yymmdd.xlsx`
 - Sheet: `ActiveEnzymes` with the following information

| rxn_id                 | enzyme_id               | direction        | kcat_values   | kegg_id               | Reactants                | Products      | EC                                               | GPR                                             | gene                                                                          | Length                       | molMass                           |
|------------------------|-------------------------|------------------|---------------|-----------------------|--------------------------|---------------|--------------------------------------------------|-------------------------------------------------|-------------------------------------------------------------------------------|------------------------------|-----------------------------------|
| str | str| Literal['f','b'] | float | str | List[str]                | List[str]     | str| str| List[str]| int| int |
|model reaction id | UniProt protein id |                  | in 1/s | reaction KEGG id | all reactants in KEGG id | all products in KEGG id | all EC numbers associated with the reaction | and/or relation between genes and reaction | list of genes coding for all peptides associated with the protein | length of protein in aa | molar mass of protein in kDa |


### 2.2 Parametrizing the TranslationalProteinSector
**Notebook**: `Scripts/i1_preprocessing/0_translational_sector_config.ipynb` (personlization needed)


**Input**:
 - Genome-scale model
 - PAM set up method (see the [PAModelpy documentation](https://pamodelpy.readthedocs.io/en/latest/)), examples are in `Modules/utils/pam_generation.py`
 - Default protein sector parameter file `Data/proteinAllocationModel_EnzymaticData_empty.xlsx` or `Results/1_preprocessing/proteinAllocationModel_EnzymaticData_<your-model>_yymmdd.xlsx` (if you already created the active enzymes sector)
 - Either one of the following:
   1. proteomics measurements of the entire proteome [g/gCDW] or [g/gPROT], an estimate of the protein coverage by the measurement, the amount of protein per gCDW, and a functional annotation of the measured peptides (e.g., with [cluster of orthologous genes](https://www.ncbi.nlm.nih.gov/research/cog/) annotations)
   2. assumption that the microorganism is similar enough to an organism which does have all from point 1. available.

**Output**: 
 - Excel file saved in `Results/1_preprocessing/proteinAllocationModel_EnzymaticData_<your-model>_yymmdd.xlsx`
 - Sheet: `Translational` with the following information for at least the relation to growth rate (`Value_for_growth` column), and possibly for the carbon source with a good amount of experimental datapoints (`Value` column):

| **Parameter**      | **Description**                                                               |
|---------------------|-------------------------------------------------------------------------------|
| `id_list`          | Identifier related to protein fraction associated with translational proteins |
| `tps_0`            | Translational protein fraction at zero growth rate. [g_protein/g_CDW]         |
| `tps_mu`           | Change in translational protein fraction per unit change of the associated reaction.  |
| `mol_mass`         | Molar mass of the translational enzymes. [kDa]                                |
| `substrate_range` | Range of values of susbstrate uptake which are associated with linear growth regime (no fermentation)|

### 2.3 Parametrizing the UnusedEnzymesSector
**Notebook**: `Scripts/i1_preprocessing/0_unused_enzyme_determination.ipynb` (personlization needed)


**Input**:
 - Default protein sector parameter file `Data/proteinAllocationModel_EnzymaticData_empty.xlsx` or `Results/1_preprocessing/proteinAllocationModel_EnzymaticData_<your-model>_yymmdd.xlsx` (if you already created the active enzymes sector)
 - Either one of the following (for the slope):
   1. maximal growth rate of a mutant obtained after adaptive laboratory evolution (ALE) in the microbes preferred carbon source
   2. assumption that the microorganism is similar enough to an organism which does have all from point 1. available. You can use the fractional increase of the growth rate after ALE to calculate the absolute maximal growth rate for your organism.
- Either one of the following (for the intercept):
   1. protein overexpression experiment in which the protein overexpression is quantified in g/g_totalprotein or g/gCDW and the growth rate for each expression strength is measured. These datapoints can be used to find the amount of unused protein at zero growth (37% for *E. coli* according to [Bruggeman et al. (2020)](https://doi.org/10.1093/femsre/fuaa034))
   2. assumption that the microorganism is similar enough to an organism which does have all from point 1. available. 

**Output**: 
 - Excel file saved in `Results/1_preprocessing/proteinAllocationModel_EnzymaticData_<your-model>_yymmdd.xlsx`
 - Sheet: `UnusedEnzyme` with the following information for at least the relation to growth rate (`Value_for_growth` column), and possibly for the carbon source with a good amount of experimental datapoints (`Value` column):

| **Parameter**     | **Description**                                                                 |
|-------------------|---------------------------------------------------------------------------------|
| `id_list`         | Identifier related to protein fraction associated with the unused enzyme sector |
| `ups_0`           | Unused enzyme fraction at zero growth rate. [g_protein/g_CDW]                   |
| `ups_mu`          | Change in unused enzyme fraction per unit change of the associated reaction.    |
| `mol_mass`        | Molar mass of unused enzymes.  [kDa]                                            |    
| `substrate_range` | Range of values of susbstrate uptake which are associated with linear growth regime (no fermentation)|

## 3. Design Choices
Now that all the parameters are parsed, there are several points to consider. In this step, it is important to keep your aim in mind.
The genetic algorithm which forms the basis of the PAMparametrizer is only informed about the diffusion limit, and not about any other 
biophysical mechanism. It is therefore important to determine to which extent the protein content affects the metabolic phenotype
of your microorganism (with the total protein content), what your starting position is (by scaling the k<sub>cat</sub> values),
which conditions you are using for parametrization, and what data you keep aside for testing the results. 

### 3.1 Total protein content
The total protein content is a key parameter as it determines how much protein the model can allocate to metabolic 
processes. It can be deduced from data or from assumptions:
- Assumption: 50% of the proteome is allocated to housekeeping proteins
- Assumption: the total protein content does not change with changing conditions
- Assumption: use the metabolic protein fraction from a closely related organism (e.g., *E. coli*: 0.258 g/gCDW)
- Quantitative proteomics data available? You can use the proteins in the model to determine the fraction of total protein mass measured which is in the model. Please correct this fraction for the amount of enzymes which is measured vs. the total amount of enzymes

For some (eukaryotic) organisms, the assumption that the total protein content does not change with changing conditions does
not hold. In these cases, you can introduce a linear relationship between the amount of protein space and any model reaction by 
adding this as a custom sector (see the [PAModelpy documentation](https://pamodelpy.readthedocs.io/en/latest/)).

### 3.2 Scaling kcat values
In some cases, the k<sub>cat</sub> values are underestimated during the initial model configuration. In this case,
it might take the genetic algorithm many iterations to increase all k<sub>cat</sub> values to prevent a too stringent protein
burden. To improve the starting position in the parametrization, you can consider to scale all k<sub>cat</sub>. You can check
which scaling factor would be most efficient for your application using the `Scripts/i1_preprocessing/1_scale_kcats_AES.py` script.

You can easily adapt the parameter file using the `increase_kcats_in_parameter_file` from PAModelpy in your PAMparametrizer setup method
 function from.

### 3.3 Which conditions?
The quality and quantity of your experimental physiology measurements determine the quality of the parametrized PAM. The PAMparametrizer
does not have any biophysical information to perform the parametrization, besides the ones you give it! If you are interested in 
a specific condition or phenotype, it makes sense to emphasize these datapoints in the parametrization by adding weights to the associated reactions.
For example in the case of *E. coli*, we are interested in the overflow metabolic phenotype, so we want to make sure that 
acetate secretion starts at the right growth rate and increases with the right slope. We can inform the PAMparametrizer about 
this by adding weights to the growth rate and the acetate exchange reaction in the HyperParams data object 
(see [PAMparamertizer setup instructions](Example.md#iii-define-the-hyperparameters)). 

In this preparatory stage, you can overlay the experimental data you have with the conditions and phenotypes you want to model.
The following questions can help you in your design choices:
 - Are there specific reactions that are representative of a specific phenotype or conditions? 
 - Are you interested in growth on different carbon sources? 
 - How much metabolic flexibility (and thus different pathways to parametrize) does your microbe have and how can you cover the entire range of pathways with your experimental measurements (e.g., different growth rates or carbon sources)?

### 3.4 Getting the validation data
Once you have decided on the conditions you want to use for parametrization, it is time to get the physiology data for 
these conditions in shape. Each carbon source gets its own ValidationData object as explained in the [PAMparamertizer setup instructions](Example.md#i-get-the-validationdata).
To easily transfer your measured fluxes to a validation data object, you need the following:
 - **Units**: mmol/gCDW/h
 - Mapping between the metabolites and the **exchange reaction ids** (such as `EX_glc__D_e` for glucose uptake)
 - **Correct sign** of measured fluxes according to the definition in the model (remember, all uptake fluxes should be negative!)
 - At least a substrate uptake rate and a growth rate per condition, preferably also O2 and CO2 exchange rates
 - All measurements for the **same strain and medium composition** (can be tricky if you are getting data from literature)

The best format to store the information is the following:

| \<Biomass reaction id\> | \<substrate uptake id\>           | \<other exchane reaction id\>     |
|-------------------------|-----------------------------------|-----------------------------------|
| ....                    | Check the direction in the model! | Check the direction in the model! |

Keep in mind that it is good practice to have separate training and validation datasets. In the end, you will need some
data to check whether the PAM you generated performs well on unseen conditions. As experimental measurements are often 
limited, it is wise to use data that the PAMparametrizer cannot use as validation data. Such as: 
 - Measurements with different cultivation strategies (autotrophy versus heterotrophy)
 - Measurements for growth on different cultivation media composition (minimal versus complex)
 - Intracellular flux measurements

## 4. Building and Running the PAMparametrizer
### 4.1 Building and running the PAMparametrizer
**Instructions**: [PAMparamertizer setup instructions](Example.md)
**Examples**: `Scripts/i2_parametrization`

**Input**:
 - Genome-scale model (see [step 1](#1-the-genome-scale-metabolic-model-sanity-check))
 - Excel file with initial sector parameters (`Results/1_preprocessing/proteinAllocationModel_EnzymaticData_<your-model>_yymmdd.xlsx`, see [step 2](#2-building-the-protein-sectors))
 - Total protein content for the PAModel (see [step 3.1](#31-total-protein-content))
 - PAM set up method (see [PAModelpy documentation](https://pamodelpy.readthedocs.io/en/latest/) and examples in `Modules/utils/pam_generation.py`)
 - Sufficient experimentally measured exchange rates (see [step 3.3](#33-which-conditions) and [step 3.4](#34-getting-the-validation-data))
 - Other design choices defined as discussed in [step 3](#3-design-choices)
 - filename extension to easily find your simulation results back, can be defined in the [HyperParams](Example.md#ii-define-the-hyperparameters) object

**Running**: pamparametrizer.run()

**Output**:
 - `Results/2_parametrization/diagnostics/pam_parametrizer_diagnostics_<your-filename-extension>.xlsx`
 - `Results/2_parametrization/progress/pam_parametrizer_progress_<your-filename-extension>.png`

### 4.2 Creating and ensemble of models
A genetic algorithm is not a deterministic method, which means the changes in k<sub>cat</sub> values are different for 
each run. Although the alternative models have a comparable phenotype, it is important to always run the framework more than once
and create what we call an 'ensemble of models'. This ensemble represents the solution space in which the turnover numbers
can deviate while still resulting in a specific metabolic phenotype. It can therefore be informative to create multiple models:
not only are you more certain of the simulated phenotype, you can also test how dependent certain modifications
to the model, strain or conditions are on the k<sub>cat</sub> values!

It is advisable to create at least 5 PAMs for your organism. The more PAMs, the more certain you can be about specific 
simulated behaviours. Nevertheless, as each parametrization can take anywhere between 5 and 30 hours on a normal desktop laptop,
5 different PAMs should sample the solution space properly without burdening your computer or HPC too much.

## 5. Analyzing the Paramaterization Results
### 5.1 The output
### 5.2 Prebuilt analyses

