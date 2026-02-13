import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
import matplotlib.colors as mcolors
from typing import List, Union, Tuple, Callable
from cobra.io import read_sbml_model

from PAModelpy import PAModel
from PAModelpy.utils import set_up_pam

from Modules.PAMparametrizer.utils.pam_generation import setup_pputida_pam, setup_cglutamicum_pam
from Modules.PAMparametrizer.utils.pam_generation import create_pamodel_from_diagnostics_file
from Modules.PAMparametrizer.utils.pamparametrizer_analysis import get_results_from_simulations

RXNS_TO_VALIDATE = {'Peripheral': ['GLCDpp', 'GAD2ktpp', 'GNK', '2DHGLCK', 'PGLCNDH', 'EX_tre_e', 'GLCt2pp'], #EX_gly_e
                    'EMP': ['HEX1', 'PGI', 'PFK', 'FBA','FBA3', 'FBP', 'TPI', 'GAPD', 'PGK', 'PGM', 'ENO', 'PYK'],
                    'ED': ['EDD', 'EDA'],
                    'PPP': ['TKT1', 'TKT2' 'G6PDH2r', 'PGL', 'GND', 'RPI', 'RPE', 'TALA'],
                    'TCA': ['PDH', 'CS', 'PC', 'OAADC', 'ACONTa', 'ACONTb', 'ICDHyr', 'SUCOAS', 'SUCDH3',
                            'SUCDi', 'FUM', 'MDH','MDH2','MDH3', 'ME2', 'ME', 'AKGDH'],
                    'Anaplerosis': ['PPC', 'PPS', 'PPCK'],
                    # 'glx': ['ICL'],
                    'Growth': ['BIOMASS_KT2440_WT3', 'BIOMASS_Ec_iML1515_core_75p37M', 'Growth'],
                    }
NEGATIVE_RXNS = ['SUCOAS', 'PGK', 'PGM', 'RPI']#reactions defined in opposite direction in the model
rxns_to_save = []
for rxns in RXNS_TO_VALIDATE.values(): rxns_to_save+=rxns

def get_simulated_fluxes_for_rxns(mfa_data:pd.Series,
                                  pam:PAModel,
                                  model_file_path: str,
                                  pam_kcat_files:List[str],
                                  setup_pam_function: Callable = set_up_pam):
    reactions_to_plot = [rxn for rxn in rxns_to_save if rxn in mfa_data]
    flux_df = pd.DataFrame(mfa_data[reactions_to_plot+['model']]).T

    m = read_sbml_model(model_file_path)
    m.optimize()

    fluxes = {'GotEnzymes': get_results_from_simulations(pam,
                                                         substrate_rates=[[mfa_data['EX_glc__D_e']]],
                                                         fluxes_to_save=reactions_to_plot,
                                                         sectors_config=False
                                                         )['fluxes'],
              'GSM': get_results_from_simulations(read_sbml_model(model_file_path),
                                                         substrate_rates=[[mfa_data['EX_glc__D_e']]],
                                                         fluxes_to_save=reactions_to_plot,
                                                         sectors_config=False
                                                         )['fluxes'],
              }

    for i, param_file in enumerate(pam_kcat_files):
        alt_pam = setup_pam_function(param_file, model = model_file_path, sensitivity=False)
        fluxes[f'Alternative {i + 1}'] = get_results_from_simulations(alt_pam,
                                                                      substrate_rates=[[mfa_data['EX_glc__D_e']]],
                                                                      fluxes_to_save=reactions_to_plot,
                                                                      sectors_config=False
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

def get_reactions2plot_pathway_mapping(flux_df):
    rxns_to_plot = {}
    for pathway, rxns in RXNS_TO_VALIDATE.items():
        rxns_in_data = [r for r in rxns if r in flux_df.columns]
        if len(rxns_in_data) > 0:
            rxns_to_plot[pathway] = rxns_in_data
    return rxns_to_plot


def get_custom_cmap_and_norm(vmin, vmax):
    # Define color segments
    colors_neg = plt.cm.Blues(np.linspace(0.4, 1, 128))
    colors_zero = np.array([[1, 1, 1, 1]])
    colors_pos = plt.cm.Reds(np.linspace(0.4, 1, 128))
    colors = np.vstack((colors_neg, colors_zero, colors_pos))
    cmap = mcolors.ListedColormap(colors, name='custom_cmap')

    # Handle 3 cases:
    if vmin < 0 and vmax > 0:
        # Case 1: both negative and positive values
        norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    elif vmax <= 0:
        # Case 2: all values are negative → use blue scale only
        cmap = plt.cm.Blues
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    else:
        # Case 3: all values are positive → use red scale only
        cmap = plt.cm.Reds
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    return cmap, norm

def plot_flux_heatmap_for_pathways(flux_df:pd.DataFrame,
                                   result_fig_path:str = None,
                                   fontsize:int=16,
                                   gs0=None,
                                   fig =None,
                                   cbar:bool =True,
                                   vrange: Tuple[int, int] = None,
                                   figsize: Tuple[int, int] = (15,10)):
    rxns_to_plot = get_reactions2plot_pathway_mapping(flux_df)
    flux_df = flux_df.rename({rxns_to_plot['Growth'][0]:'Growth rate'}, axis=1)

    if vrange is None:
        vmax = np.ceil(flux_df.max().max())
        vmin = np.floor(flux_df.min().min())
    else:
        (vmin,vmax) = vrange


    # Create heatmap
    if fig is None: fig= plt.figure(figsize=figsize)
    if gs0 is None: gs0 = gridspec.GridSpec(1,1)[0]

    plt.rcParams.update({'font.size': fontsize})

    combined_cmap, norm = get_custom_cmap_and_norm(vmin=vmin, vmax=vmax)

    #plot growth rate with different colorbar
    if 'Growth' in rxns_to_plot:
        growth_df = flux_df[['Growth rate']]
        flux_df = flux_df.drop('Growth rate', axis =1)

        gs = gridspec.GridSpecFromSubplotSpec(
            1, 6,
            subplot_spec=gs0,
            width_ratios=[5, 0.4, 0.2,0.3,0.8, 0.3],  # main heatmap, cbar, growth heatmap, cbar
            wspace=0.05
        )
        if cbar:
            ax_main = fig.add_subplot(gs[0, 0])
            cax_main = fig.add_subplot(gs[0, 3])
            ax_growth = fig.add_subplot(gs[0, 1])
            cax_growth = fig.add_subplot(gs[0, 5])

            sns.heatmap(growth_df, ax=ax_growth, cmap="Greys", cbar=cbar, cbar_ax=cax_growth,
                        yticklabels=False, xticklabels=True,
                        vmin=0.2, vmax = 0.6)
                        # vmin = 0.3, vmax = 0.6)
                        # vmin=round(growth_df.min(),1),
                        # vmax=round(growth_df.max(),1))

            # ax_growth.set_title(r"Growth rate", fontsize = fontsize, fontweight = 'bold')
            # ax_growth.set_ylabel(r"Growth rate [$\text{h}^{-1}$]", fontsize=fontsize)
            ax_growth.tick_params(labelsize=fontsize)  # Set tick label font size
            ax_growth.set_xticklabels(ax_growth.get_xticklabels(), rotation=45, ha='right')
            # ax_growth.tick_params(axis='x', which = 'both', bottom = False, labelbottom = False, top=False)
            ax_growth.set_ylabel("")

            sns.heatmap(flux_df, annot=False, cmap=combined_cmap, norm=norm, fmt=".2f", ax=ax_main,
                                  vmax=vmax, vmin=vmin, cbar=cbar, cbar_ax = cax_main
                                  )
            # ax_main.tick_params(axis='y', labelrotation=90)
            # ax_main.set_yticklabels(ax_main.get_yticklabels(), rotation=90, ha='right')

            cbar_growth = ax_growth.collections[0].colorbar
            cbar_growth.ax.tick_params(labelsize=fontsize)  # Set tick label font size
            cbar_growth.set_label('Growth rate [1/h]', fontsize=fontsize)

        else:
            # ax.tick_params(axis='y', labelrotation=90)
            ax_main = fig.add_subplot(gs[0, :5])
            ax_growth = fig.add_subplot(gs[0, 5])

            sns.heatmap(growth_df, ax=ax_growth, cmap="Greys", cbar=cbar,
                        xticklabels=True, vmin=0.2, vmax=0.6)
            # vmin=round(growth_df.min(), 1),
            # vmax=round(growth_df.max(), 1))

            ax_growth.set_ylabel("")  # no duplicate
            ax_growth.set_xticklabels(ax_growth.get_xticklabels(), rotation=45, ha='right')
            ax_growth.tick_params(axis='x', labelsize=fontsize - 1)
            ax_growth.set_yticks([])
            ax_growth.set_yticklabels([])
            # ax_last.tick_params(axis='y', left=False)

            sns.heatmap(flux_df, annot=False, cmap=combined_cmap, norm=norm, fmt=".2f", ax=ax_main,
                                  vmax=vmax, vmin=vmin, cbar=cbar,
                                  )

        for a in [ax_main, ax_growth]:
            for side in ['right', 'left', 'top', 'bottom']:
                a.spines[side].set_visible(False)


    else:
        ax_main = fig.subplots(gs0)
        # ax.tick_params(axis='y', labelrotation=90)
        ax_main.tick_params(axis='both', labelsize=fontsize)
        sns.heatmap(flux_df, annot=False, cmap=combined_cmap, norm=norm, fmt=".2f", ax=ax_main,
                vmax = vmax, vmin = vmin, cbar = cbar
                )

    # draw line after experimental data
    ax_main.axhline(y=1, color='black', linewidth=2)
    ax_main.set_xticklabels(ax_main.get_xticklabels(), rotation=45, ha='right')
    ax_main.set_ylabel("")
    ax_main.tick_params(axis='both', labelsize=fontsize-1)

    # Access and modify the colorbar
    if cbar:
        cbar = ax_main.collections[0].colorbar
        cbar.ax.tick_params(labelsize=fontsize)  # Set tick label font size
        cbar.set_label(r'Flux [mmol/$\text{g}_{\text{CDW}}$/h]', fontsize=fontsize)

    # add pathway annotation
    group_labels = []
    group_positions = []
    start = 0
    for group_name, columns in rxns_to_plot.items():
        if group_name == 'Growth': continue
        length = len(columns)
        center = start + length / 2
        group_labels.append(group_name)
        group_positions.append(center)
        start += length
        # drawing line between the groups
        if start < flux_df.shape[1]:  # Avoid drawing a line at the far right
            ax_main.axvline(x=start, color='black', linewidth=2)

    # Add group labels as a second x-axis (on top)
    ax_top = ax_main.twiny()
    # Match the ticks and limits
    ax_top.set_xlim(ax_main.get_xlim())
    ax_top.set_xticks(group_positions)
    ax_top.set_xticklabels(group_labels, fontsize=fontsize, fontweight='bold')
    ax_top.tick_params(axis='x', bottom=False, top=True, labelbottom=False, labeltop=True)
    ax_top.spines['bottom'].set_visible(False)


    if result_fig_path is not None:
        print(f'Saving figure to {result_fig_path}')
        plt.tight_layout()
        # plt.show()
        plt.savefig(result_fig_path)
    else:
        return ax_main

def main_iJN1463(gs = None, fig=None,cbar=True, vrange= None, fontsize = 11):
    NUM_MODELS = 5
    PAM_KCAT_FILES_IJN = [os.path.join('Results', '3_analysis', 'parameter_files',
                                       f'proteinAllocationModel_EnzymaticData_iJN1463_{file_nmbr}.xlsx') for file_nmbr in
                          range(1, NUM_MODELS + 1)]
    PPUTIDA_PHENOTYPE_FILE_PATH = os.path.join('Data', 'Pputida_phenotypes', 'pputida_phenotypes.xlsx')

    mfa_data = pd.read_excel(PPUTIDA_PHENOTYPE_FILE_PATH,
                             sheet_name='fluxomics_glucose').iloc[0]  # first row has most relevant data
    mfa_data['model'] = 'Nikel, et al (2015)'
    pam = setup_pputida_pam(sensitivity=False)

    flux_df = get_simulated_fluxes_for_rxns(mfa_data,
                                            pam,
                                            os.path.join('Models', 'iJN1463.xml'),
                                            PAM_KCAT_FILES_IJN,
                                            setup_pam_function=setup_pputida_pam)
    fig_out = None
    if gs is None: fig_out = os.path.join('Results', '3_analysis', 'Metabolic_pathways_pputida.png')
    ax = plot_flux_heatmap_for_pathways(flux_df,
                                   result_fig_path=fig_out,
                                   gs0=gs,
                                   fig=fig,
                                   cbar=cbar,
                                   vrange = vrange,
                                        figsize=(25,10),
                                        fontsize=fontsize
                                   )
    return ax


def main_iML1515(gs=None, fig = None, fig_out=None, cbar=True, vrange = None, num_models: int = 10):
    NUM_MODELS = num_models
    PAM_KCAT_FILES_IML = [os.path.join('Results', '3_analysis', 'parameter_files',
                                       f'proteinAllocationModel_EnzymaticData_iML1515_{file_nmbr}.xlsx') for file_nmbr in
                          range(1, NUM_MODELS + 1)]
    PAM_KCAT_FILES_IML.append(os.path.join('Results', '3_analysis', 'parameter_files',
                                       'proteinAllocationModel_EnzymaticData_iML1515_csources.xlsx'))

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
                                            os.path.join('Models', 'iML1515.xml'),
                                            PAM_KCAT_FILES_IML)
    if gs is None and fig_out is None:
        fig_out = os.path.join('Results', '3_analysis', 'Metabolic_pathways_ecoli.png')
    ax = plot_flux_heatmap_for_pathways(flux_df,
                                   result_fig_path=fig_out,
                                   gs0=gs,
                                   fig=fig,
                                   cbar=cbar,
                                   vrange = vrange,
                                   )
    return ax


def main_iCGB21FR(gs=None, fig = None, cbar=True, vrange = (0, 8), fontsize = 11):
    NUM_MODELS = 5
    PAM_KCAT_FILES_ICGB = [os.path.join('Results', '3_analysis', 'parameter_files',
                                       f'proteinAllocationModel_EnzymaticData_iCGB21FR_{file_nmbr}.xlsx') for file_nmbr in
                          range(1, NUM_MODELS + 1)]
    cglutamicum_phenotype_file_path = os.path.join('Data', 'Cglutamicum_phenotypes', 'cglutamicum_phenotypes.xlsx')

    mfa_data = pd.read_excel(cglutamicum_phenotype_file_path,
                             sheet_name='fluxomics_glucose').iloc[1]  # second row has the flux data in the right units
    mfa_data['model'] = 'Peifer, et al (2012)'
    pam = setup_cglutamicum_pam(sensitivity=False)

    flux_df = get_simulated_fluxes_for_rxns(mfa_data,
                                            pam,
                                            os.path.join('Models', 'iCGB21FR_annotated_copyable.xml'),
                                            PAM_KCAT_FILES_ICGB,
                                            setup_pam_function=setup_cglutamicum_pam)
    fig_out = None
    if gs is None:fig_out = os.path.join('Results', '3_analysis', 'Metabolic_pathways_cglutamicum.png')
    ax = plot_flux_heatmap_for_pathways(flux_df,
                                   result_fig_path=fig_out,
                                   gs0=gs,
                                   fig=fig,
                                   vrange = vrange,
                                   cbar = cbar,fontsize=fontsize)
    return ax

if __name__ == '__main__':
    main_iML1515()
    # main_iJN1463()
    # main_iCGB21FR()



