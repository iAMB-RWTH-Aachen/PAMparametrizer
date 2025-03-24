import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from Modules.utils.pam_generation import setup_pputida_pam
from Modules.utils.pam_generation import create_pamodel_from_diagnostics_file
from Modules.utils.pamparametrizer_analysis import get_results_from_simulations

NUM_MODELS = 5

PAM_KCAT_FILES_IJN = [os.path.join('Results', '2_parametrization', 'diagnostics',
                               f'pam_parametrizer_diagnostics_iJN1463_{file_nmbr}.xlsx') for file_nmbr in
                  range(1, NUM_MODELS + 1)]

PPUTIDA_PHENOTYPE_FILE_PATH = os.path.join('Data', 'Pputida_phenotypes', 'pputida_phenotypes.xlsx')
RXNS_TO_VALIDATE = {'Peripheral': ['GLCDpp', 'GAD2ktpp','GNK', '2DHGLCK', 'PGLCNDH'],
                    'EMP': ['HEX1','PGI', 'FBA', 'FBP', 'GAPD', 'PGK', 'PGM', 'ENO', 'PYK'],
                    'ED': ['EDD','EDA'],
                    'PPP': ['TKT1', 'G6PDH2r', 'PGL','GND', 'RPI', 'RPE', 'TALA'],
                    'TCA': ['PDH', 'PC','OAADC','ACONTa','ACONTb','ICDHyr','SUCOAS','SUCDi','FUM','MDH','ME2'],
                    'Growth': ['BIOMASS_KT2440_WT3']
                    # 'Anaplerosis':['PPC'],
                    # 'glx shunt':['ICL']
                    }
rxns_to_save = []
for rxns in RXNS_TO_VALIDATE.values(): rxns_to_save+=rxns

def get_simulated_fluxes_for_rxns():
    mfa_data = pd.read_excel(PPUTIDA_PHENOTYPE_FILE_PATH,
                             sheet_name='fluxomics_glucose').iloc[0]  # first row has most relevant data
    flux_df = pd.DataFrame(mfa_data[rxns_to_save]).T
    flux_df['model'] = 'Nikel, et al (2015)'
    pam = setup_pputida_pam(sensitivity=False)

    fluxes = {'GotEnzymes': get_results_from_simulations(pam,
                                                         substrate_rates=[[mfa_data['EX_glc__D_e']]],
                                                         fluxes_to_save=rxns_to_save,
                                                         transl_sector_config=False
                                                         )['fluxes']
              }

    for i, alternative_file in enumerate(PAM_KCAT_FILES_IJN):
        alt_pam = create_pamodel_from_diagnostics_file(
            alternative_file,
            pam.copy(copy_with_pickle=True),
            other_enzyme_id_pattern=r'E[0-9][0-9]*|Enzyme_PP_\d+'  # find the pputida enzyme ids created from locus tags

        )
        fluxes[f'Alternative {i + 1}'] = get_results_from_simulations(alt_pam,
                                                                      substrate_rates=[[mfa_data['EX_glc__D_e']]],
                                                                      fluxes_to_save=rxns_to_save,
                                                                      transl_sector_config=False
                                                                      )['fluxes']
    for model, df in fluxes.items():
        df['model'] = model
        flux_df = pd.concat([flux_df, df])

    flux_df.to_excel('Results/3_analysis/mfa_experiment_pputida.xlsx')
    return flux_df

def plot_flux_heatmap_for_pathways(flux_df:pd.DataFrame):
    flux_df = flux_df.drop(['substrate_id', 'substrate'], axis=1).set_index('model').astype(float)

    # Create heatmap
    fig, ax = plt.subplots(figsize=(20, 20))
    sns.heatmap(flux_df, annot=True, cmap="coolwarm", fmt=".2f",
                ax=ax, cbar_kws={'label': r'flux [mmol/$\text{g}_{\text{CDW}}$/h]'}
                )
    # draw line after experimental data
    ax.axhline(y=1, color='black', linewidth=2)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

    # add pathway annotation
    # Prepare group labels and positions
    group_labels = []
    group_positions = []

    start = 0
    for group_name, columns in RXNS_TO_VALIDATE.items():
        length = len(columns)
        center = start + length / 2
        group_labels.append(group_name)
        group_positions.append(center)
        start += length
        # drawing line between the groups
        if start < flux_df.shape[1]:  # Avoid drawing a line at the far right
            ax.axvline(x=start, color='black', linewidth=2)

    # Add group labels as a second x-axis (on top)
    ax_top = ax.twiny()

    # Match the ticks and limits
    ax_top.set_xlim(ax.get_xlim())
    ax_top.set_xticks(group_positions)
    ax_top.set_xticklabels(group_labels, fontsize=10, fontweight='bold')
    ax_top.tick_params(axis='x', bottom=False, top=True, labelbottom=False, labeltop=True)
    ax_top.spines['bottom'].set_visible(False)

    plt.title("Flux comparison heatmap")
    plt.xlabel("Reaction")
    plt.ylabel("Model")
    plt.tight_layout()
    plt.savefig(os.path.join('Results', '3_analysis', 'Metabolic_pathways_pputida.png'))


if __name__ == '__main__':
    flux_df = get_simulated_fluxes_for_rxns()
    plot_flux_heatmap_for_pathways(flux_df)




