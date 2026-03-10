import os
from typing import List, Dict
import pandas as pd
from matplotlib import gridspec
import seaborn as sns
import matplotlib.pyplot as plt
from cobra.io import read_sbml_model

from Modules.PAMparametrizer.utils.pam_generation import setup_cglutamicum_pam, setup_pputida_pam
from PAModelpy import PAModel
from Figures.Scripts.Figure1_iml1515_kcat_analysis import recreate_progress_plot
from Scripts.i2_parametrization.pam_parametrizer_iCGB21FR import set_up_validation_data as refdata_setup_icgb21fr
from Scripts.i2_parametrization.pam_parametrizer_iJN1463 import set_up_validation_data as refdata_setup_ijn1463


from Scripts.i3_analysis.metabolic_flux_distribution_vs_exp import main_iJN1463 as plot_intracell_flux_distribution_ijn
from Scripts.i3_analysis.metabolic_flux_distribution_vs_exp import main_iCGB21FR as plot_intracell_flux_distribution_icgb

from Modules.PAMparametrizer.utils.pam_generation import create_pamodel_from_diagnostics_file

NUM_MODELS = 5
PAM_KCAT_FILES_ICG = [os.path.join('Results', '2_parametrization', 'diagnostics',
                               f'pam_parametrizer_diagnostics_iCGB21FR_{file_nmbr}.xlsx') for file_nmbr in
                  range(1, NUM_MODELS + 1)]

PAM_KCAT_FILES_IJN = [os.path.join('Results', '2_parametrization', 'diagnostics',
                               f'pam_parametrizer_diagnostics_iJN1463_{file_nmbr}.xlsx') for file_nmbr in
                  range(1, NUM_MODELS + 1)]
FONTSIZE=11

labels = [f'Alternative {i}' for i in range(1, NUM_MODELS+1)]
rxns2label = {'EX_o2_e': r'O$_{2}$ uptake'+'\n'+r'[mmol/g$_{\text{CDW}}$/h]',
              'EX_co2_e': r'CO$_{2}$ evolution'+'\n'+r'[mmol/g$_{\text{CDW}}$/h]',
              'EX_glc__D_e': 'glc uptake'+'\n'+r'[mmol/g$_{\text{CDW}}$/h]',
              'BIOMASS_KT2440_WT3':'Growth [1/h]', 'Growth':'Growth [1/h]'}

def plot_simulations_vs_experiments(pamodel: 'PAModel',
                                    diagnostic_files: List[str],
                                    other_enzyme_id_pattern: str,
                                    gem:'Model',
                                    exp_data: pd.DataFrame,
                                    to_plot: List[str],
                                    cmap: Dict[str, 'Color'],
                                    fig: plt.Figure, gs: 'GridSpec',
                                    sub_uptake: str = 'EX_glc__D_e',
                                    enzyme_sector_update: bool = False):
    models = {'GEM': gem, 'After preprocessing':pamodel}
    for i, file in enumerate(diagnostic_files):
        alt_pam = create_pamodel_from_diagnostics_file(file,
                                                       pamodel.copy(copy_with_pickle =True),
                                                       enzyme_sector_update=enzyme_sector_update,
                                                       other_enzyme_id_pattern = other_enzyme_id_pattern
                                                       )
        models[f'Alternative {i+1}'] = alt_pam
    axs = [fig.add_subplot(gs[j]) for j in range(len(to_plot))]

    for j, rxn in enumerate(to_plot):
        axs[j].scatter([abs(r) for r in exp_data[sub_uptake]], [abs(f) for f in exp_data[rxn]], color='black')
        axs[j].grid(visible=True, alpha=0.2, linewidth=0.7)
        axs[j].set_axisbelow(True)

    for model_id, model in models.items():
        sub_rates = []
        flux = []
        for rate in exp_data[sub_uptake]:
            if isinstance(model, PAModel):
                model.change_reaction_bounds(sub_uptake, rate, 0)
            else:
                model.reactions.get_by_id(sub_uptake).lower_bound = rate
            sol = model.optimize()
            if model.solver.status == 'optimal':
                sub_rates+= [abs(rate)]
                flux += [sol.fluxes]
        for j,rxn in enumerate(to_plot):
            kwargs = {'color':cmap[model_id]}
            if model_id == 'GEM': kwargs = {'linestyle':'--', 'color': 'black'}
            axs[j].plot(sub_rates, [abs(f[rxn]) for f in flux], **kwargs, label = model_id)

    return axs


def main():
    model_colors = sns.color_palette("viridis", n_colors=NUM_MODELS)
    other_colors = {'GotEnzymes': 'grey', 'After preprocessing': 'black', 'iCGB21FR': 'purple','iJN1463': 'purple', 'GEM': 'grey'}
    cmap = {
        **{l: c for l, c in
           zip([f'Alternative {i + 1}' for i in range(NUM_MODELS)], model_colors)},
        **other_colors}

    pam_icgb = setup_cglutamicum_pam(os.path.join('Results', '2_parametrization',
                                                  'proteinAllocationModel_iCGB21FR_EnzymaticData_multi.xlsx'),
                                     sensitivity=False
                                     )
    pam_ijn = setup_pputida_pam(os.path.join('Results', '2_parametrization',
                                                  'proteinAllocationModel_iJN1463_EnzymaticData_multi.xlsx'),
                                     sensitivity=False
                                     )
    exp_data_icgb = refdata_setup_icgb21fr(pam_info_file=os.path.join('Results', '2_parametrization',
                                                  'proteinAllocationModel_iCGB21FR_EnzymaticData_multi.xlsx'),
                                           csources = ['Glucose'])[0].valid_data.sort_values('EX_glc__D_e')
    exp_data_ijn = refdata_setup_ijn1463(pam_info_file=os.path.join('Results', '2_parametrization',
                                                  'proteinAllocationModel_iJN1463_EnzymaticData_multi.xlsx'),
                                           csources = ['Glucose'])[0].valid_data.sort_values('EX_glc__D_e')
    fig = plt.figure(figsize=(21/2.56, 30/2.56))
    gs_main = gridspec.GridSpec(2, 1, height_ratios=[1,1], hspace=0.45)
    gs_main_top = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[1,1], hspace=0.6,subplot_spec=gs_main[0])
    gs_main_bottom = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[4,3], hspace=0.6, subplot_spec=gs_main[1])

    gs_cgb_fluxes = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs_main_top[0],
                                                    wspace=0.5)

    gs_cgb_heatmap = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_main_top[1],width_ratios=[1,20])[1]
    # gs_cgb = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs_main[1])
    gs_ijn_fluxes_legend =gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_main_bottom[-1], width_ratios=[2,2])
    gs_ijn_heatmap = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_main_bottom[0],width_ratios=[1,20])[1]

    for title, gs in zip(
            [r'$\mathbf{\it{Corynebacterium\, glutamicum}}$ ATCC 13032', r'$\mathbf{\it{Pseudomonas\, putida}}$ KT2440'],
            [gs_main[0],gs_main[1]]
    ):
        ax = fig.add_subplot(gs)
        ax.axis('off')
        ax.set_title(title, pad = 35, fontsize = FONTSIZE+4)

    i=0
    axs = []
    for pamodel, gs,gem_file, kcat_file_list, regex, rxns_to_plot, exp_data, update_enz_sector in zip(
            [pam_icgb, pam_ijn],
            [gs_cgb_fluxes, gs_ijn_fluxes_legend],
            [os.path.join('Models', file) for file in ['iCGB21FR_annotated_copyable.xml', 'iJN1463.xml']],
            [PAM_KCAT_FILES_ICG, PAM_KCAT_FILES_IJN],
            [r'Enzyme_cg[0-9]+',r'Enzyme_*|Enzyme_PP_[0-9]+'],
            [['Growth', 'EX_co2_e', 'EX_o2_e'],['BIOMASS_KT2440_WT3']],
            [exp_data_icgb, exp_data_ijn],
            [False, False]
    ):
        gem = read_sbml_model(gem_file)
        for rxn in [gem.reactions.EX_o2_e, gem.reactions.EX_co2_e]:#in icgb21FR oxygen and co2 uptake are bounded
            rxn.bounds = (-1000,1000)
        ax = plot_simulations_vs_experiments(pamodel = pamodel,
                                             gem=gem,
                                             gs=gs, fig = fig,
                                             diagnostic_files=kcat_file_list,
                                             other_enzyme_id_pattern = regex,
                                             to_plot=rxns_to_plot,
                                             cmap = cmap,
                                             exp_data=exp_data,
                                             enzyme_sector_update=update_enz_sector)
        for i,a in enumerate(ax):
            a.set_ylabel(rxns2label[rxns_to_plot[i]])
        axs+=ax
        i+=1

    plot_intracell_flux_distribution_icgb(gs = gs_cgb_heatmap, fig = fig, fontsize = FONTSIZE, vrange=(-6,13))
    plot_intracell_flux_distribution_ijn(gs = gs_ijn_heatmap, fig = fig, fontsize = FONTSIZE, vrange=(-6,13), cbar = False)

    legend_ax = fig.add_subplot(gs_ijn_fluxes_legend[1])
    legend_ax.axis("off")  # Hide axes
    h, l = axs[0].get_legend_handles_labels()

    legend_ax.legend(h, l, loc="center",
                     fontsize=FONTSIZE,
                     ncol=2, frameon=False)
    #add annotation
    annotations = ["","","A", "B", "C", "F", "D", "","", "","","","","E"]

    for ax, label in zip(fig.axes, annotations):
        ax.annotate(label, xy=(-0.1, 1), xycoords="axes fraction",
                    fontsize=FONTSIZE, fontweight='bold',
                    xytext=(-5, 5), textcoords="offset points",
                    ha="right", va="bottom")

    # Row 0: Centered across all 4 axes
    # fig.text(-0.1, 0.75, 'Corynebacterium glutamicum', ha='center', va='center', fontsize=FONTSIZE, weight = 'bold', rotation = 90)
    # fig.text(-0.1, 0.25, 'Pseudomonas putida', ha='center', va='center', fontsize=FONTSIZE, weight = 'bold', rotation = 90)
    for x,y in zip([(fig.axes[2].get_position().x0+fig.axes[4].get_position().x1)/2,
                    (fig.axes[5].get_position().x0+fig.axes[5].get_position().x1)/2],
                   [fig.axes[2].get_position().y0+0.005,fig.axes[5].get_position().y0]
                   ):
        fig.text(x,y-0.04, r'Glucose uptake rate [mmol/$\text{g}_\text{CDW}$/h]', ha='center', va='center', fontsize=FONTSIZE)


    # plt.tight_layout()
    fig.tight_layout()
    fig.savefig(os.path.join('Figures', 'Figure5_cglutamicum_pputida.png'))

if __name__ == '__main__':
    main()