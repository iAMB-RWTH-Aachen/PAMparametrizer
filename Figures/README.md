# Figures associated with publication related to this repository

This repository is related to the following publication:....
For reproducibility, the code associated with the figures in this publication can be found in this directory. 
The code for some supplementary figures are associated with standard analysis scripts in `Scripts/i3_analysis`.
The following table gives an overview which figures and tables are generated with which code:


| Figure or Table                 | Description                                                                                                           | Where can I find the code?                                                                              |
|---------------------------------|-----------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|
| **Main Figures**                |
| Figure 1 and Table S3           | Analysis of parametrization results of the *E. coli* PAM                                                              | `Figure/Scripts/Figure1_iml1515_kcat_analysis.py`                                                       |
| Figure 2, Table S4 and Table S5 | Errors and sensitivities of the *E. coli* PAMs                                                                        | `Figure/Scripts/Figure2_sensitivity_error.py`                                                           |
| Figure 3                        | Parametrization of the *C. glutanicum* and  *P. putida* PAMs                                                          | `Figure/Scripts/Figure3_alternative_models.py`                                                          |
| Figure 4                        | Flux distribution the *C. glutanicum* and  *P. putida* PAMs and effects of data reduction                             | `Figure/Scripts/Figure4_data_reduction.py`                                                              |
| **Supplementary Figures**       |
| Figure S3 and Table S1          | Parametrization of the translational protein sector of the *E. coli* PAM                                              | `Scripts/i1_preprocessing/0_translational_sector_config.ipynb`                                          |
| Figure S4                       | Error progression of the PAMparametrizer and genetic algorithm                                                        | `Figure/Scripts/SuppFigures.py`                                                                         |
| Figure S5 and Table S3          | Cumulative flux distribution of the *E. coli* PAMs                                                                    | `Figure/Scripts/SuppFigures.py`                                                                         |
| Figure S6                       | Distribution of protein concentrations among COG for *E. coli* PAMs                                                   | `Scripts/i3_analysis/pam_parametrizer_validate_proteomics.ipynb`                                        |
| Figure S7                       | Kernel densitity plot of protein concentrations for *E. coli* PAMs                                                    | `Scripts/i3_analysis/pam_parametrizer_validate_proteomics.ipynb`                                        |
| Figure S8                       | Simulated intracellular fluxes for growth on glucose at several dilution rates for *E. coli* PAMs                     | `Scripts/i3_analysis/PAMparametrizer_iML1515_validation.ipynb`                                          |
| Figure S9                       | Simulated intracellular fluxes for growth on glucose for *E. coli* PAMs (heatmap)                                     | `Scripts/i3_analysis/PAMparametrizer_iML1515_validation.ipynb`                                          |
| Figure S10 and Table S4         | Boxplots of difference between simulated and measured fluxes forgrowth on different carbon sources for *E. coli* PAMs | `Scripts/i3_analysis/PAMparametrizer_iML1515_validation_csources.ipynb`                                 |
| Figure S11                      | Simulated intracellular fluxes for growth on different carbon sources for *E. coli* PAMs                              | `Scripts/i3_analysis/PAMparametrizer_iML1515_validation_csources.ipynb`                                 |
| **Supplementary Tables**        |
| Table S1                        | Parametrization of the protein sector  of the *E. coli*                                                               | `Scripts/i1_preprocessing/0_translational_sector_config.ipynb`/   `0_unused_enzyme_determination.ipynb` |
| Table S7                        | Parametrization of the protein sector for *C. glutanicum*                                                             | `Scripts/i1_preprocessing/translational_sector_config_iCGB21FR.ipynb` |
