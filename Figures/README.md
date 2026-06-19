# Figures associated with publication related to this repository

This repository is related to the following publication:....
For reproducibility, the code associated with the figures in this publication can be found in this directory. 
The code for some supplementary figures are associated with standard analysis scripts in `Scripts/i3_analysis`.
The following table gives an overview which figures and tables are generated with which code:


| Figure or Table                 | Description                                                                                                                                   | Where can I find the code?                                                                              |
|---------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|
| **Main Figures**                |
| Figure 2                        | Simulated exchange rates for 10 alternative *E. coli* PAMs as a function of substrate uptake rate                                             | `Figures/Scripts/Figure2_alternative_pams.py`                                                           |
| Figure 3 and Table S4           | Analysis of parametrization results of the *E. coli* PAM                                                                                      | `Figures/Scripts/Figure3_kcat_protein_analysis.py`                                                      |
| Figure 4, Table S4 and Table S5 | Errors and sensitivities of the *E. coli* PAMs                                                                                                | `Figures/Scripts/Figure4_sensitivity_error.py`                                                          |
| Figure 5 and Table S10          | Parametrization of the *C. glutanicum* and  *P. putida* PAMs                                                                                  | `Figures/Scripts/Figure4_alternative_models.py`                                                         |
| Figure 6                        | The effects of data reduction on parametrization performance on the *E. coli* PAM                                                             | `Figures/Scripts/Figure6_data_reduction.py`                                                             |
| **Supplementary Figures**       |
| Figure S3 and Table S1          | Parametrization of the translational protein sector of the *E. coli* PAM                                                                      | `Scripts/i1_preprocessing/0_translational_sector_config.ipynb`                                          |
| Figure S4                       | Error progression of the PAMparametrizer and genetic algorithm                                                                                | `Figures/Scripts/SuppFigures.py`                                                                        |
| Figure S5 and Table S4          | Effect of initial parameters on the prediction accuracy of PAMs before and after parametrization                                              | `Figures/Scripts/SuppFig_differen_model_analysis.py`                                                    |
| Figure S6 and Table S5          | Cumulative flux distribution of the *E. coli* PAMs                                                                                            | `Figures/Scripts/SuppFigures.py`                                                                        |
| Figure S7                       | Distribution of protein concentrations among COG for *E. coli* PAMs                                                                           | `Scripts/i3_analysis/pam_parametrizer_validate_proteomics.ipynb`                                        |
| Figure S8                       | Overview of coeffivients of variation of the K<sub>cat<\sub> values, flux rates, and protein concentrations in the alternative *E. coli* PAMs | `Figures/Scripts/flux_and_kcat_variability`                                                             |
| Figure S9                       | Simulated intracellular fluxes for growth on glucose at several dilution rates for *E. coli* PAMs                                             | `Scripts/i3_analysis/PAMparametrizer_iML1515_validation.ipynb`                                          |
| Figure S10                      | Simulated intracellular fluxes for growth on glucose for *E. coli* PAMs (heatmap)                                                             | `Scripts/i3_analysis/PAMparametrizer_iML1515_validation.ipynb`                                          |
| Figure S11 and Table S6 and S9  | Boxplots of difference between simulated and measured fluxes for growth on different carbon sources for *E. coli* PAMs                        | `Scripts/i3_analysis/PAMparametrizer_iML1515_validation_csources.ipynb`                                 |
| Figure S12                      | Simulated intracellular fluxes for growth on different carbon sources for *E. coli* PAMs                                                      | `Scripts/i3_analysis/PAMparametrizer_iML1515_validation_csources.ipynb`                                 |
| **Supplementary Tables**        |
| Table S1                        | Parametrization of the protein sector  of the *E. coli*                                                                                       | `Scripts/i1_preprocessing/0_translational_sector_config.ipynb`/   `0_unused_enzyme_determination.ipynb` |
| Table S11                       | Parametrization of the protein sector for *C. glutanicum*                                                                                     | `Scripts/i1_preprocessing/translational_sector_config_iCGB21FR.ipynb`                                   |
