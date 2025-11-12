import os

import pandas as pd
from matplotlib import gridspec
import seaborn as sns
import matplotlib.pyplot as plt
from cobra.io import read_sbml_model

from Modules.utils.pam_generation import setup_cglutanicum_pam, setup_pputida_pam
from PAModelpy import PAModel
from Figures.Scripts.Figure1_iml1515_kcat_analysis import recreate_progress_plot
from Scripts.i2_parametrization.pam_parametrizer_iCGB21FR import set_up_validation_data as refdata_setup_icgb21fr
from Scripts.i2_parametrization.pam_parametrizer_iJN1463 import set_up_validation_data as refdata_setup_ijn1463


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
rxns2label = {'EX_o2_e': 'O2 uptake\n[mmol/gCDW/h]', 'EX_co2_e': 'CO2 excretion\n[mmol/gCDW/h]', 'EX_glc__D_e': 'glc uptake\n[mmol/gCDW/h]', 'Growth': 'Growth',
              'BIOMASS_KT2440_WT3':'Growth [1/h]'}

def plot_simulations_vs_experiments(pamodel: 'PAModel', diagnostic_files, gem,
                                    exp_data: pd.DataFrame, to_plot,cmap,
                                    fig, gs,sub_uptake: str = 'EX_glc__D_e'):
    models = {'GEM': gem, 'After preprocessing':pamodel}
    i = 5
    for file in diagnostic_files:
        models[f'Alternative {i}'] = create_pamodel_from_diagnostics_file(file,
                                                                          pamodel.copy(copy_with_pickle =True),
                                                                          enzyme_sector_update=False)
    axs = [fig.add_subplot(gs[j]) for j in range(len(to_plot))]
    print(models)

    for j, rxn in enumerate(to_plot):
        axs[j].scatter([abs(r) for r in exp_data[sub_uptake]], [abs(f) for f in exp_data[rxn]], color='black')
        axs[j].grid()

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
            axs[j].plot(sub_rates, [abs(f[rxn]) for f in flux], **kwargs)
    return axs





def main():
    # fig, axs = plt.subplots(figsize=(20,20))

    # gs = gridspec.GridSpec(nrows=2,ncols=1,height_ratios=[3,2], figure=fig)
    # gs_top = gridspec.GridSpecFromSubplotSpec(ncols = 4, nrows=2, subplot_spec=gs[0])
    # gs_bottom = gridspec.GridSpecFromSubplotSpec(ncols = 2, nrows=1, subplot_spec=gs[1])

    model_colors = sns.color_palette("viridis", n_colors=NUM_MODELS)
    other_colors = {'GotEnzymes': 'grey', 'After preprocessing': 'black', 'iCGB21FR': 'purple','iJN1463': 'purple', 'GEM': 'grey'}
    cmap = {
        **{l: c for l, c in
           zip([f'Alternative {i + 1}' for i in range(NUM_MODELS)], model_colors)},
        **other_colors}

    pam_icgb = setup_cglutanicum_pam(os.path.join('Results', '2_parametrization',
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
    gs_main = gridspec.GridSpec(4, 1, height_ratios=[4,4,4,3], hspace=0.8)

    gs_cgb_fluxes = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs_main[0],
                                                    wspace=0.4)

    gs_cgb_heatmap = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_main[1],width_ratios=[1,10])[1]
    # gs_cgb = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs_main[1])
    gs_ijn_fluxes_legend =gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_main[-1], width_ratios=[2,2])
    gs_ijn_heatmap = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_main[2],width_ratios=[1,10])[1]

    i=0
    for pamodel, gs,gem_file, kcat_file_list, rxns_to_plot, exp_data in zip(
            [pam_icgb, pam_ijn],
            [gs_cgb_fluxes, gs_ijn_fluxes_legend],
            [os.path.join('Models', file) for file in ['iCGB21FR_annotated_copyable.xml', 'iJN1463.xml']],
            [PAM_KCAT_FILES_ICG, PAM_KCAT_FILES_IJN],
            [['Growth', 'EX_co2_e', 'EX_o2_e'],['BIOMASS_KT2440_WT3']],
            [exp_data_icgb, exp_data_ijn]
    ):
        gem = read_sbml_model(gem_file)
        for rxn in [gem.reactions.EX_o2_e, gem.reactions.EX_co2_e]:#in icgb21FR oxygen and co2 uptake are bounded
            rxn.bounds = (-1000,1000)
        # ax = [fig.add_subplot(gs[j]) for j in range(len(rxns_to_plot))]
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
        axs= plot_simulations_vs_experiments(pamodel = pamodel,
                                        gem=gem,
                                        gs=gs, fig = fig,
                                        diagnostic_files=kcat_file_list,
                                        to_plot=rxns_to_plot,
                                        cmap = cmap,
                                        exp_data=exp_data)
        for i,a in enumerate(axs):
            a.set_ylabel(rxns2label[rxns_to_plot[i]])

        i+=1

    plot_intracell_flux_distribution_icgb(gs = gs_cgb_heatmap, fig = fig, fontsize = FONTSIZE, vrange=(-6,13))
    plot_intracell_flux_distribution_ijn(gs = gs_ijn_heatmap, fig = fig, fontsize = FONTSIZE, vrange=(-6,13), cbar = False)

    legend_ax = fig.add_subplot(gs_ijn_fluxes_legend[-1])
    legend_ax.axis("off")  # Hide axes
    h, l = fig.axes[1].get_legend_handles_labels()

    legend_ax.legend(h, l, loc="center",
                     fontsize=FONTSIZE,
                     ncol=2, frameon=False)
    #add annotation
    annotations = ["A", "B", "C", "F", "D", "","", "","","","","E"]

    for ax, label in zip(fig.axes, annotations):
        ax.annotate(label, xy=(-0.1, 1), xycoords="axes fraction",
                    fontsize=FONTSIZE, fontweight='bold',
                    xytext=(-5, 5), textcoords="offset points",
                    ha="right", va="bottom")

    # Row 0: Centered across all 4 axes
    # fig.text(-0.1, 0.75, 'Corynebacterium glutanicum', ha='center', va='center', fontsize=FONTSIZE, weight = 'bold', rotation = 90)
    # fig.text(-0.1, 0.25, 'Pseudomonas putida', ha='center', va='center', fontsize=FONTSIZE, weight = 'bold', rotation = 90)
    for x,y in zip([(fig.axes[0].get_position().x0+fig.axes[2].get_position().x1)/2,
                    (fig.axes[-2].get_position().x0+fig.axes[-2].get_position().x1)/2],
                   [fig.axes[0].get_position().y0,fig.axes[-1].get_position().y0]
                   ):
        fig.text(x,y-0.05, r'Glucose uptake rate [mmol/$\text{g}_\text{CDW}$/h]', ha='center', va='center', fontsize=FONTSIZE)


    # plt.tight_layout()
    fig.tight_layout()
    fig.savefig(os.path.join('Figures', 'Figure3_cglutanicum_pputida.png'))

if __name__ == '__main__':
    main()