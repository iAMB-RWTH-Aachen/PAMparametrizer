import os

import pandas as pd
from matplotlib import gridspec
import seaborn as sns
import matplotlib.pyplot as plt
from PAModelpy import PAModel
from Figures.Scripts.Figure1_iml1515_kcat_analysis import recreate_progress_plot
from Scripts.i2_parametrization.pam_parametrizer_iCGB21FR import set_up_pamparametrizer as pamparam_setup_icgb21fr
from Scripts.i2_parametrization.pam_parametrizer_iJN1463 import set_up_pamparametrizer as pamparam_setup_ijn1463


from Scripts.i3_analysis.metabolic_flux_distribution_vs_exp import main_iJN1463 as plot_intracell_flux_distribution_ijn
from Scripts.i3_analysis.metabolic_flux_distribution_vs_exp import main_iCGB21FR as plot_intracell_flux_distribution_icgb

from Modules.utils.pam_generation import create_pamodel_from_diagnostics_file

NUM_MODELS = 5
PAM_KCAT_FILES_ICG = [os.path.join('Results', '2_parametrization', 'diagnostics',
                               f'pam_parametrizer_diagnostics_iCGB21FR_{file_nmbr}.xlsx') for file_nmbr in
                  range(1, NUM_MODELS + 1)]

PAM_KCAT_FILES_IJN = [os.path.join('Results', '2_parametrization', 'diagnostics',
                               f'pam_parametrizer_diagnostics_iJN1463_{file_nmbr}.xlsx') for file_nmbr in
                  range(1, NUM_MODELS -1)]
FONTSIZE=11

labels = [f'Alternative {i}' for i in range(1, NUM_MODELS+1)]
rxns2label = {'EX_o2_e': 'O2 uptake', 'EX_co2_e': 'CO2 excretion', 'EX_glc__D_e': 'glc uptake', 'Growth': 'Growth',
              'BIOMASS_KT2440_WT3':'Growth'}

def plot_simulations_vs_experiments(pamodel: 'PAModel', diagnostic_files, gem,
                                    exp_data: pd.DataFrame, to_plot,
                                    fig, gs,sub_uptake: str = 'EX_glc__D_e'):
    models = {'GEM': gem, 'After preprocessing':pamodel}
    i = 1
    for file in diagnostic_files:
        models[f'Alternative {i}'] = create_pamodel_from_diagnostics_file(file,
                                                                          pamodel.copy(copy_with_pickle =True),
                                                                          enzyme_sector_update=False)
    ax_flux = fig.add_subplot(gs)

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
                sub_rates+= [rate]
                flux += [sol['fluxes']]
        for rxn in to_plot:
            if (('biomass' in rxn.lower()) or ('growth' in rxn.lower())) and len(to_plot)>1:pass
            ax_flux.plot(sub_rates, [abs(f[rxn]) for f in flux])





def main():
    # fig, axs = plt.subplots(figsize=(20,20))

    # gs = gridspec.GridSpec(nrows=2,ncols=1,height_ratios=[3,2], figure=fig)
    # gs_top = gridspec.GridSpecFromSubplotSpec(ncols = 4, nrows=2, subplot_spec=gs[0])
    # gs_bottom = gridspec.GridSpecFromSubplotSpec(ncols = 2, nrows=1, subplot_spec=gs[1])

    model_colors = sns.color_palette("viridis", n_colors=NUM_MODELS)
    other_colors = {'GotEnzymes': 'grey', 'After preprocessing': 'black', 'iCGB21FR': 'purple','iJN1463': 'purple'}
    cmap = {
        **{l: c for l, c in
           zip([f'Alternative {i + 1}' for i in range(NUM_MODELS)], model_colors)},
        **other_colors}

    fig = plt.figure(figsize=(21/2.56, 30/2.56))
    gs_main = gridspec.GridSpec(4, 1, height_ratios=[4,4,4,3], hspace=0.5)

    gs_cgb_fluxes = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs_main[0],
                                                    wspace=0.4)
    # gs_cgb = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs_main[1])
    gs_ijn_fluxes_legend =gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_main[-1], width_ratios=[2,1])
    i=0
    axs = []
    for pamparamsetup, gs,gem_file, kcat_file_list, kwargs, rxns_to_plot in zip(
            [pamparam_setup_icgb21fr, pamparam_setup_ijn1463],
            [gs_cgb_fluxes, gs_ijn_fluxes_legend],
            [os.path.join('Models', file) for file in ['iCGB21FR_annotated_copyable.xml', 'iJN1463.xml']],
            [PAM_KCAT_FILES_ICG, PAM_KCAT_FILES_IJN],
            [{'max_substrate_uptake_rate':-0.1,
              'c_sources': ['Glucose', 'Fructose', 'Succinate','Gluconate', 'Acetate'],
                'kcat_increase_factor': 6
              },
             {'max_substrate_uptake_rate':-0.1,
              'min_substrate_uptake_rate':-15,
              'c_sources': ['Glycerol', 'Glucose','Octanoate','m-Xylene','Succinate', 'Benzoate', 'Fructose'],
                'kcat_increase_factor': 5
              }
             ],
            [['Growth', 'EX_co2_e', 'EX_o2_e'],['BIOMASS_KT2440_WT3']]
    ):
        ax = [fig.add_subplot(gs[j]) for j in range(len(rxns_to_plot))]
        # ax = axs[i]
        # recreate_progress_plot(kcat_file_list,
        #                         labels, fig, ax,
        #                         legend = False,
        #                         fontsize=FONTSIZE,
        #                         pamparam_setup=pamparamsetup,
        #                         pamparam_kwargs = kwargs,
        #                         rxns_to_plot = rxns_to_plot,
        #                        gem_file = gem_file,
        #                        cmap = cmap,
        #                         other_measurements = False,
        #                        enzyme_sector_update = False)
        axs.append(ax)
        # for i,a in enumerate(ax[:-1]):
        #     a.set_ylabel(rxns2label[rxns_to_plot[i]])

        i+=1


    plot_intracell_flux_distribution_icgb(gs = gs_main[1], fig = fig, fontsize = FONTSIZE, vrange=(-6,13))
    plot_intracell_flux_distribution_ijn(gs = gs_main[2], fig = fig, fontsize = FONTSIZE, vrange=(-6,13), cbar = False)

    legend_ax = fig.add_subplot(gs_ijn_fluxes_legend[-1])
    legend_ax.axis("off")  # Hide axes
    h, l = fig.axes[0].get_legend_handles_labels()

    legend_ax.legend(h, l, loc="center",
                     fontsize=FONTSIZE,
                     ncol=1, frameon=False)
    #add annotation
    annotations = ["A", "B", "C", "D", "F", "G","E", "","","","","H"]

    for ax, label in zip(fig.axes, annotations):
        ax.annotate(label, xy=(-0.1, 1), xycoords="axes fraction",
                    fontsize=FONTSIZE, fontweight='bold',
                    xytext=(-5, 5), textcoords="offset points",
                    ha="right", va="bottom")

    # Row 0: Centered across all 4 axes
    # fig.text(-0.1, 0.75, 'Corynebacterium glutanicum', ha='center', va='center', fontsize=FONTSIZE, weight = 'bold', rotation = 90)
    # fig.text(-0.1, 0.25, 'Pseudomonas putida', ha='center', va='center', fontsize=FONTSIZE, weight = 'bold', rotation = 90)
    for x,y in zip([(fig.axes[0].get_position().x0+fig.axes[2].get_position().x1)/2, (fig.axes[5].get_position().x0+fig.axes[6].get_position().x1)/2],[0.75,0.66]):
        fig.text(x,y, r'Glucose uptake rate [mmol/$\text{g}_\text{CDW}$/h]', ha='center', va='center', fontsize=FONTSIZE)


    # plt.tight_layout()
    fig.tight_layout()
    fig.savefig(os.path.join('Figures', 'Figure3_cglutanicum_pputida.png'))

if __name__ == '__main__':
    main()