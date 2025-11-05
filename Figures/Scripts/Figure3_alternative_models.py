import os
from matplotlib import gridspec
import seaborn as sns
import matplotlib.pyplot as plt

from Figures.Scripts.Figure1_iml1515_kcat_analysis import recreate_progress_plot
from Scripts.i2_parametrization.pam_parametrizer_iCGB21FR import set_up_pamparametrizer as pamparam_setup_icgb21fr
from Scripts.i2_parametrization.pam_parametrizer_iJN1463 import set_up_pamparametrizer as pamparam_setup_ijn1463


from Scripts.i3_analysis.metabolic_flux_distribution_vs_exp import main_iJN1463 as plot_intracell_flux_distribution_ijn
from Scripts.i3_analysis.metabolic_flux_distribution_vs_exp import main_iCGB21FR as plot_intracell_flux_distribution_icgb


NUM_MODELS = 5
PAM_KCAT_FILES_ICG = [os.path.join('Results', '2_parametrization', 'diagnostics',
                               f'pam_parametrizer_diagnostics_iCGB21FR_{file_nmbr}.xlsx') for file_nmbr in
                  range(1, NUM_MODELS + 1)]

PAM_KCAT_FILES_IJN = [os.path.join('Results', '2_parametrization', 'diagnostics',
                               f'pam_parametrizer_diagnostics_iJN1463_{file_nmbr}.xlsx') for file_nmbr in
                  range(1, NUM_MODELS -1)]
FONTSIZE=11

labels = [f'Alternative {i}' for i in range(1, NUM_MODELS+1)]

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

    # create a pretty figure
    fig = plt.figure(figsize=(30/2.56, 20/2.56))

    # Outer GridSpec: 2 rows (80% heatmaps, 20% colorbar)
    gs_main = gridspec.GridSpec(3, 1, height_ratios=[4,4,0.2])
    gs_cgb = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_main[0],
                                                  wspace=0.2, width_ratios=[1, 1]
                                                  )

    gs_cgb_fluxes = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs_cgb[0, 0],
                                                    wspace=0.4,
                                                    hspace=0.4)
    gs_ijn = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_main[1],
                                                  wspace=0.2, width_ratios=[1, 4]
                                                  )

    gs_ijn_fluxes = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_ijn[0, 0],
                                                    hspace=0.2)

    i=0
    for pamparamsetup, gs,gem_file, kcat_file_list, kwargs, rxns_to_plot in zip(
            [pamparam_setup_icgb21fr, pamparam_setup_ijn1463],
            [gs_cgb_fluxes, gs_ijn_fluxes],
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
        ax = [fig.add_subplot(gs[j]) for j in range(len(rxns_to_plot)+1)]
        # ax = axs[i]
        recreate_progress_plot(kcat_file_list,
                                labels, fig, ax,
                                legend = False,
                                fontsize=FONTSIZE,
                                pamparam_setup=pamparamsetup,
                                pamparam_kwargs = kwargs,
                                rxns_to_plot = rxns_to_plot,
                               gem_file = gem_file,
                               cmap = cmap,
                                other_measurements = True,
                               enzyme_sector_update = False)
        i+=1


    plot_intracell_flux_distribution_icgb(gs = gs_cgb[1], fig = fig, fontsize = FONTSIZE, vrange=(-6,13))
    plot_intracell_flux_distribution_ijn(gs = gs_ijn[1], fig = fig, fontsize = FONTSIZE, vrange=(-6,13), cbar = False)

    legend_ax = fig.add_subplot(gs_main[2])
    legend_ax.axis("off")  # Hide axes
    h, l = fig.axes[0].get_legend_handles_labels()

    legend_ax.legend(h, l, loc="center",
                     fontsize=FONTSIZE,
                     ncol=5, frameon=False)
    #add annotation
    annotations = ["A", "B", "C", "D", "E", "F", "G", "H"]

    for ax, label in zip(fig.axes, annotations):
        ax.annotate(label, xy=(0, 1), xycoords="axes fraction",
                    fontsize=FONTSIZE, fontweight='bold',
                    xytext=(-5, 5), textcoords="offset points",
                    ha="right", va="bottom")

    # Row 0: Centered across all 4 axes
    # left = fig.axes[1, 0].get_position().x0
    # right = fig.axes[1, 2].get_position().x1
    # mid = (left + right) / 2
    # fig.text(0.5, 0.95, 'Corynebacterium glutanicum', ha='center', va='center', fontsize=FONTSIZE, weight = 'bold')
    # fig.text(mid, 0.52, r'Glucose uptake rate [mmol/$\text{g}_\text{CDW}$/h]', ha='center', va='center', fontsize=FONTSIZE)
    #
    # # Row 1: Centered between axs[1,1] and axs[1,2]
    # # We'll find the x-position of those two axes and take their midpoint
    # left = axs[1, 0].get_position().x0
    # right = axs[1, 1].get_position().x1
    # mid = (left + right) / 2
    #
    # fig.text(mid, 0.48, 'Pseudomonas putida', ha='center', va='center', fontsize=FONTSIZE, weight = 'bold')
    # fig.text((axs[1, 0].get_position().x0+axs[1, 0].get_position().x1)/2,
    #          0.05, r'Glucose uptake rate [mmol/$\text{g}_\text{CDW}$/h]', ha='center', va='center', fontsize=FONTSIZE)


    plt.tight_layout()
    # fig.tight_layout()
    fig.savefig(os.path.join('Figures', 'Figure3_cglutanicum_pputida.png'))

if __name__ == '__main__':
    main()