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
- \>60% of all reactions should be annotated with gene-reaction associations

The following is not required, but makes your life easier:
- Reactions annotated with KEGG ids
- Metabolites annotated with KEGG ids
- Genes annotated with KEGG ids

Also, make sure you have the following information
- Regular expression to recognize gene-ids (e.g., `r'\b([b|s]\d{4})\b'` for the b1234 locus tags in *E. coli*)
- Potential exceptions in identifier/mapping format (an example is `s0001` which is used for unknown gene)

## 2. Building the protein sectors
You can find scripts to help you build the protein sectors in `Scripts/i1_preprocessing`. The scripts can be run in random
order, although we suggest starting with the ActiveEnzymesSector.

### 2.1 Initiating the ActiveEnzymesSector
Notebook: `Scripts/i1_preprocessing/0_parse_kcat_values_GotEnzymes.ipynb` (personlization needed)
Input:
 - genome-scale model
 - default protein sector parameter file `Data/proteinAllocationModel_EnzymaticData_empty.xlsx`
 - **UniProt** mapping between genes and proteins (instructions in notebook)
 - All k<sub>cat<\sub>s for the microorganism from **GotEnzymes** (instructions in notebook)
Output: 
 - Excel file saved in `Results/1_preprocessing/proteinAllocationModel_EnzymaticData_<your-model>_yymmdd.xlsx`
 - Sheet: `ActiveEnzymes` with the following information

| rxn_id                 | enzyme_id               | direction         | kcat_values   | kegg_id               | Reactants                     |Products| EC                                               | GPR                                             | gene                                                                          | Length                       | molMass                           |
|------------------------|-------------------------|-------------------|---------------|-----------------------|-------------------------------|-----|--------------------------------------------------|-------------------------------------------------|-------------------------------------------------------------------------------|------------------------------|-----------------------------------|
| str: model reaction id | str: UniProt protein id | Literal['f', 'b'] | float: in 1/s | str: reaction KEGG id | List[str]: all KEGG compounds |List[str]: all KEGG compounds| str: all EC numbers associated with the reaction | str: and/or relation between genes and reaction | List[str] : list of genes coding for all peptides associated with the protein | int: length of protein in aa | int: molar mass of protein in kDa |


### 2.2 Parametrizing the TranslationalProteinSector
Notebook: `Scripts/i1_preprocessing/0_translational_sector_config.ipynb` (personlization needed)
Input:
 - genome-scale model
 - pam set up method (see the [PAModelpy](https://github.com/iAMB-RWTH-Aachen/PAModelpy) documentation)
 - default protein sector parameter file `Data/proteinAllocationModel_EnzymaticData_empty.xlsx` or `Results/1_preprocessing/proteinAllocationModel_EnzymaticData_<your-model>_yymmdd.xlsx` (if you already created the active enzymes sector)
 - either one of the following:
   1. proteomics measurements of the entire proteome [g/gCDW] or [g/gPROT], an estimate of the protein coverage by the measurement, to amount of protein per gCDW
   2. assumption that the microorganism is similar enough to an organism which does have all from point 1. available.
Output: 
 - Excel file saved in `Results/1_preprocessing/proteinAllocationModel_EnzymaticData_<your-model>_yymmdd.xlsx`
 - Sheet: `Translational` with the following information

### 2.3 Parametrizing the UnusedEnzymesSector

## 3. Design Choices
### 3.1 Total protein content
### 3.2 Scaling kcat values

### 3.3 Which conditions?
### 3.3 Getting the validation data

## 4. Building and Running the PAMparametrizer
### 4.1 Building the PAMparametrizer
### 4.2 Creating and ensemble of models

## 5. Analyzing the Paramaterization Results
### 5.1 The output
### 5.2 Prebuild analyses

