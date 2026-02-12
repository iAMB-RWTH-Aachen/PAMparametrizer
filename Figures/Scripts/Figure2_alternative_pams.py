import matplotlib.pyplot as plt
from cobra.io.sbml import read_sbml_model
from matplotlib.colors import to_hex
from matplotlib import gridspec
import seaborn as sns
import numpy as np
import pandas as pd
import os
from typing import List, Dict, Tuple, Callable
from PAModelpy.utils import set_up_pam

from Scripts.i2_parametrization.pam_parametrizer_iML1515 import set_up_pamparametrizer

from Modules.PAMparametrizer.utils.pam_generation import (get_rxn2kcat_protein2gene_dict,
                                          _extract_reaction_id_from_catalytic_reaction_id,
                                          _get_rxn2kcat_as_series,
                                          create_pamodel_from_diagnostics_file
                                         )
from Modules.PAMparametrizer.utils.pamparametrizer_analysis import (calculate_kcat_differences,
                                                    plot_histogram_logspace,
                                                    select_clustered_rows_by_variation,
                                                    get_clusters_from_clustermap,
                                                    set_up_pam_parametrizer_and_get_substrate_uptake_rates
                                                   )

from Modules.PAMparametrizer.utils.pamparametrizer_visualization import plot_valid_data, plot_simulation, plot_flux_vs_experiment
from Modules.PAMparametrizer.utils.pamparametrizer_setup import set_up_sector_config_from_diagnostic_file

def recreate_progress_plot(best_indiv_files:list[str],
                           labels:list[str], fig,
                           axs,
                           legend:bool = True,
                           fontsize=20,
                           pamparam_setup: Callable= None,
                           pamparam_kwargs: dict = {'max_substrate_uptake_rate':-0.1},
                           rxns_to_plot: List[str] = None,
                           substrate_uptake_id:str = 'EX_glc__D_e',
                           gem_file: str = os.path.join('Models', 'iML1515.xml'),
                           other_measurements: bool = False,
                           cmap:dict = None,
                           enzyme_sector_update:bool = True
                           ):
    rxn_mapper = {'BIOMASS_Ec_iML1515_core_75p37M': 'Growth rate [1/h]',
                  'EX_co2_e': r'CO$_2$ excretion',
                  'EX_o2_e': 'rO$_2$ uptake',
                  'EX_ac_e': 'Acetate excretion'}
    j=0

    if pamparam_setup is None:
        parametrizer, substrate_rates = set_up_ecoli_pam_parametrizer_and_get_substrate_uptake_rates()
    else:
        parametrizer, substrate_rates = set_up_pam_parametrizer_and_get_substrate_uptake_rates(pamparam_setup,
                                                                                               pamparam_kwargs)
    if rxns_to_plot is not None:
        parametrizer.validation_data.get_by_id(substrate_uptake_id)._reactions_to_plot = rxns_to_plot
    else:
        rxns_to_plot = parametrizer.validation_data.get_by_id(substrate_uptake_id)._reactions_to_plot

    substrate_rates = sorted(substrate_rates)
    gem = read_sbml_model(gem_file)
    fig, axs = plot_valid_data(parametrizer,axs, fig, fontsize=fontsize)
    print('Run reference simulations')

    gem_fluxes = []
    for rate in substrate_rates:
        gem.reactions.EX_glc__D_e.lower_bound = rate
        sol = gem.optimize()
        gem_fluxes.append(sol.fluxes)
    pam = parametrizer.pamodel_no_sensitivity.copy(copy_with_pickle=True)


    fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                                   substrate_rates = substrate_rates,
                                                                   sensitivity = False)
    fig, axs = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates],
                               parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
                               iteration=0, color='black',label = 'After \npreprocessing')

    fig, axs = plot_simulation(fig, axs, gem_fluxes, [abs(rate) for rate in substrate_rates],
                               parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
                               iteration=0, color='black', label='iML1515', plotting_kwargs = {'linestyle':'--'})

    for file, label in zip(best_indiv_files, labels):
        j +=1
        print('\nAlternative ', label, ' from file ', file)
        parametrizer.pamodel_no_sensitivity = create_pamodel_from_diagnostics_file(file,
                                                                    pam.copy(copy_with_pickle = True),
                                                                    enzyme_sector_update = enzyme_sector_update)
        parametrizer.pamodel_no_sensitivity.change_reaction_bounds('EX_glc__D_e', -5,0)
        if enzyme_sector_update:
            parametrizer.validation_data.EX_glc__D_e.sector_configs = set_up_sector_config_from_diagnostic_file(file)

        fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                                       substrate_rates=substrate_rates,
                                                                       sensitivity=False)
        fig, axs, color = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates],
                                   parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
                                   iteration=j + 1, max_iteration=len(best_indiv_files), label = label,
                                          color = cmap[label],return_color=True)
        if other_measurements:
            print('plotting other carbon sources')
            plot_flux_vs_experiment(axs[len(rxns_to_plot)], parametrizer,
                                    color, fontsize)

    for rxn, ax in zip(rxns_to_plot,axs):
        ax.set_ylabel(rxn_mapper[rxn])
    lines, labels = fig.axes[1].get_legend_handles_labels()

    if legend:
        fig.legend(lines, labels, loc='upper left', bbox_to_anchor=(0.1, 0.99), frameon=False,
                   fontsize=fontsize - 5)
    # fig.tight_layout()
    return fig, axs

def set_up_ecoli_pam_parametrizer_and_get_substrate_uptake_rates() -> Tuple:
    kwargs = {'min_substrate_uptake_rate':-12,
              'max_substrate_uptake_rate': -0.1,
              'kcat_increase_factor': 7}
    return set_up_pam_parametrizer_and_get_substrate_uptake_rates(set_up_pamparametrizer,
                                                           kwargs)

def main():
    NUM_ALT_MODELS = 10
    FONTSIZE = 14
    diagnostic_files = [os.path.join('Results', '2_parametrization', 'diagnostics',
                                     f'pam_parametrizer_diagnostics_{file_nmbr}.xlsx') for file_nmbr in
                        range(1, NUM_ALT_MODELS + 1)]

    model_colors = sns.color_palette("viridis", n_colors=NUM_ALT_MODELS)
    other_colors = {'GotEnzymes': 'grey', 'After \npreprocessing': 'black', 'iML1515': 'purple'}
    cmap = {
        **{l: c for l, c in
           zip([f'Alternative {i + 1}' for i in range(NUM_ALT_MODELS)], model_colors)},
        **other_colors}

    #create a pretty figure
    fig = plt.figure(figsize=(21/2.58, 17/2.58))

    gs_main = gridspec.GridSpec(1, 2,wspace=0.6,
                                width_ratios=[10,1.2])
    gs_inner = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs_main[0],
                                                    wspace = 0.4,
                                                    hspace=0.4
                                                )

    axs1 =[None,None,None,None]
    for i in reversed(range(4)):  # Assuming line_axs[0] is a row of axes
        row, col = (i, 0)  # Determine row/column position in the 2x2 grid
        kwargs={}
        if i>1:
            row, col = (i-2, 1)
        if i==0 or i==1:
            kwargs={'sharex':axs1[i+1]}
        axs1[i]= fig.add_subplot(gs_inner[row, col], **kwargs)
    line, line_axs = recreate_progress_plot(diagnostic_files,
                                           labels=[f'Alternative {i}' for i in range(1, NUM_ALT_MODELS + 1)],
                                            fig=fig, axs=axs1, legend=False, fontsize=FONTSIZE,
                                            cmap = cmap)
    #share x and y axis labels for progress plot
    ax_group = fig.add_subplot(gs_main[0])
    ax_group.set_xticks([])
    ax_group.set_yticks([])
    ax_group.set_frame_on(False)
    ax_group.set_ylabel(r"Flux [mmol/$\text{g}_{\text{CDW}}$/h]",
                        labelpad=20, fontsize=FONTSIZE)
    ax_group.yaxis.set_label_coords(-0.15, 0.5)
    ax_group.set_xlabel(r"Glucose uptake [mmol/$\text{g}_{\text{CDW}}$/h]",
                        labelpad=20, fontsize=FONTSIZE)
    # create a legend
    legend_ax = fig.add_subplot(gs_main[1])
    # for ax in [line_axs[0], ax2]:
    legend_ax.axis("off")  # Hide axes
    # line_axs[0].plot([],[],label='After preprocessing', color='black', linewidth=5)#dummy line for complete legend
    line_axs[0].plot([],[],label='iML1515', color='black', linestyle ='--',linewidth=5)#dummy line for complete legend

    handles, labels = line_axs[0].get_legend_handles_labels()
    h, l = [],[]
    for handle, label in zip(handles, labels):
        if not label in l:
            l.append(label)
            h.append(handle)

    legend_ax.legend(h, l, loc="center",
                     fontsize=FONTSIZE,
                     ncol=1, frameon=False)
    # ax2.legend(handles, labels, loc="lower center", ncol=round(len(labels)/2), frameon=False)
    for ax in line_axs:
        ax.grid(visible=True, alpha=0.2, linewidth=0.7)
        ax.set_axisbelow(True)

    #Add alphabet annotations
    #first 4 axes are added in reverse, and must ignore additional axis for x/y label config
    annotations = ["D", "C", "B", "A",""]
    fontsize = FONTSIZE  # Adjust as needed

    for ax, label in zip(fig.axes, annotations):
        ax.annotate(label, xy=(0, 1), xycoords="axes fraction",
                    fontsize=fontsize, fontweight='bold',
                    xytext=(-5, 5), textcoords="offset points",
                    ha="right", va="bottom")
    fig.tight_layout()
    fig.savefig(os.path.join('Figures', 'Figure2_parametrization_results_flux_analysis.png'))

if __name__ == '__main__':
    main()