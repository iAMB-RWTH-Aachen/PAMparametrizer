import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Union

from PAModelpy import PAModel
from PAModelpy.utils import set_up_pam

from Modules.utils.pam_generation import setup_pputida_pam, setup_cglutanicum_pam
from Modules.utils.pam_generation import create_pamodel_from_diagnostics_file
from Modules.utils.pamparametrizer_analysis import get_results_from_simulations

RXNS_TO_VALIDATE = {'Peripheral': ['GLCDpp', 'GAD2ktpp', 'GNK', '2DHGLCK', 'PGLCNDH'],
                    'EMP': ['HEX1', 'PGI', 'PFK', 'FBA', 'FBP', 'TPI', 'GAPD', 'PGK', 'PGM', 'ENO', 'PYK'],
                    'ED': ['EDD', 'EDA'],
                    'PPP': ['TKT1', 'TKT2' 'G6PDH2r', 'PGL', 'GND', 'RPI', 'RPE', 'TALA'],
                    'TCA': ['PDH', 'CS', 'PC', 'OAADC', 'ACONTa', 'ACONTb', 'ICDHyr', 'SUCOAS', 'SUCDH3',
                            'SUCDi', 'FUM', 'MDH', 'ME2', 'ME'],
                    'Anaplerosis': ['PPC', 'PPS', 'PPCK'],
                    'glx shunt': ['ICL'],
                    'Growth': ['BIOMASS_KT2440_WT3', 'BIOMASS_Ec_iML1515_core_75p37M', 'Growth'],
                    }
NEGATIVE_RXNS = ['SUCOAS', 'PGK', 'PGM']#reactions defined in opposite direction in the model
rxns_to_save = []
for rxns in RXNS_TO_VALIDATE.values(): rxns_to_save+=rxns

def get_simulated_fluxes_for_rxns(mfa_data:pd.Series,
                                  pam:PAModel,
                                  pam_kcat_files:List[str]):
    reactions_to_plot = [rxn for rxn in rxns_to_save if rxn in mfa_data]
    flux_df = pd.DataFrame(mfa_data[reactions_to_plot]).T

    fluxes = {'GotEnzymes': get_results_from_simulations(pam,
                                                         substrate_rates=[[mfa_data['EX_glc__D_e']]],
                                                         fluxes_to_save=reactions_to_plot,
                                                         transl_sector_config=False
                                                         )['fluxes']
              }

    for i, alternative_file in enumerate(pam_kcat_files):
        alt_pam = create_pamodel_from_diagnostics_file(
            alternative_file,
            pam.copy(copy_with_pickle=True),
            # find the pputida cglutanicum enzyme ids created from locus tags or reaction_ids
            other_enzyme_id_pattern=r'E[0-9][0-9]*|Enzyme_PP_\d+|Enzyme_cg\d+|Enzyme_\w+_(\d+.)(\d+.)(\d+.)(\d+])|Enzyme_\w+'

        )
        fluxes[f'Alternative {i + 1}'] = get_results_from_simulations(alt_pam,
                                                                      substrate_rates=[[mfa_data['EX_glc__D_e']]],
                                                                      fluxes_to_save=reactions_to_plot,
                                                                      transl_sector_config=False
                                                                      )['fluxes']
    for model, df in fluxes.items():
        #make sure the sign of the flux rates align
        for neg_rxn in NEGATIVE_RXNS:
            if neg_rxn in df.columns: df[neg_rxn] = -df[neg_rxn]
        df['model'] = model
        flux_df = pd.concat([flux_df, df])

    flux_df = flux_df.drop(['substrate_id', 'substrate'], axis=1).set_index('model').astype(float)
    flux_df = flux_df.dropna(axis=1, how='any')

    return flux_df

def plot_flux_heatmap_for_pathways(flux_df:pd.DataFrame,
                                   result_fig_path:str,
                                   vmax:Union[int,float]=None,
                                   fontsize:int=16):
    rxns_to_plot = {}
    for pathway, rxns in RXNS_TO_VALIDATE.items():
        rxns_in_data = [r for r in rxns if r in flux_df.columns]
        if len(rxns_in_data)>0:
            rxns_to_plot[pathway] = rxns_in_data
    if vmax is None:
        vmax = np.ceil(flux_df.max().max())
        print(vmax)

    # Create heatmap
    fig, ax = plt.subplots(figsize=(20, 12))
    heatmap = sns.heatmap(flux_df, annot=False, cmap="coolwarm", fmt=".2f",
                ax=ax, cbar_kws={'label': r'flux [mmol/$\text{g}_{\text{CDW}}$/h]'},
                vmax = vmax
                )
    # draw line after experimental data
    ax.axhline(y=1, color='black', linewidth=2)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_yticklabels(ax.get_xticklabels(), rotation=90, ha='right')
    ax.tick_params(axis='both', labelsize=fontsize)

    # Access and modify the colorbar
    cbar = heatmap.collections[0].colorbar
    cbar.ax.tick_params(labelsize=fontsize)  # Set tick label font size

    # add pathway annotation
    group_labels = []
    group_positions = []
    start = 0
    for group_name, columns in rxns_to_plot.items():
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

    # plt.title("Flux comparison heatmap")
    ax.set_xlabel("Reaction", fontsize=fontsize)
    plt.ylabel("Model", fontsize=fontsize)
    plt.tight_layout()
    plt.savefig(result_fig_path)

def main_iJN1463():
    NUM_MODELS = 5
    PAM_KCAT_FILES_IJN = [os.path.join('Results', '2_parametrization', 'diagnostics',
                                       f'pam_parametrizer_diagnostics_iJN1463_{file_nmbr}.xlsx') for file_nmbr in
                          range(1, NUM_MODELS + 1)]
    PPUTIDA_PHENOTYPE_FILE_PATH = os.path.join('Data', 'Pputida_phenotypes', 'pputida_phenotypes.xlsx')

    mfa_data = pd.read_excel(PPUTIDA_PHENOTYPE_FILE_PATH,
                             sheet_name='fluxomics_glucose').iloc[0]  # first row has most relevant data
    mfa_data['model'] = 'Nikel, et al (2015)'
    pam = setup_pputida_pam(sensitivity=False)

    flux_df = get_simulated_fluxes_for_rxns(mfa_data,
                                            pam,
                                            PAM_KCAT_FILES_IJN)

    plot_flux_heatmap_for_pathways(flux_df,
                                   os.path.join('Results', '3_analysis', 'Metabolic_pathways_pputida.png'))


def main_iML1515():
    NUM_MODELS = 8
    PAM_KCAT_FILES_IML = [os.path.join('Results', '2_parametrization', 'diagnostics',
                                       f'pam_parametrizer_diagnostics_{file_nmbr}.xlsx') for file_nmbr in
                          range(1, NUM_MODELS + 1)]

    mfa_data = pd.read_excel(os.path.join('Data', 'Ecoli_phenotypes', 'fluxomics_datasets_gerosa.xlsx'))
    mfa_data_glc = mfa_data[mfa_data.condition == 'Glucose'][['reaction', 'measured']].set_index('reaction').squeeze()
    mfa_data_glc['model'] = 'Gerosa, et al (2015)'
    pam = set_up_pam(os.path.join('Results',
                                  '2_parametrization',
                                  'proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx'
                                  )
                     ,sensitivity=False
                     )

    flux_df = get_simulated_fluxes_for_rxns(mfa_data_glc,
                                            pam,
                                            PAM_KCAT_FILES_IML)
    plot_flux_heatmap_for_pathways(flux_df,
                                   os.path.join('Results', '3_analysis', 'Metabolic_pathways_ecoli.png'))


def main_iCGB21FR():
    NUM_MODELS = 5
    PAM_KCAT_FILES_ICGB = [os.path.join('Results', '2_parametrization', 'diagnostics',
                                       f'pam_parametrizer_diagnostics_iCGB21FR_{file_nmbr}.xlsx') for file_nmbr in
                          range(1, NUM_MODELS + 1)]
    cglutnicum_phenotype_file_path = os.path.join('Data', 'Cglutanicum_phenotypes', 'cglutanicum_phenotypes.xlsx')

    mfa_data = pd.read_excel(cglutnicum_phenotype_file_path,
                             sheet_name='fluxomics_glucose').iloc[1]  # second row has the flux data in the right units
    mfa_data['model'] = 'Peifer, et al (2012)'
    pam = setup_cglutanicum_pam(sensitivity=False)

    flux_df = get_simulated_fluxes_for_rxns(mfa_data,
                                            pam,
                                            PAM_KCAT_FILES_ICGB)
    plot_flux_heatmap_for_pathways(flux_df,
                                   os.path.join('Results', '3_analysis', 'Metabolic_pathways_cglutanicum.png'),
                                   vmax = 8)

if __name__ == '__main__':
    main_iML1515()
    main_iJN1463()
    main_iCGB21FR()



