import os
from matplotlib import gridspec
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import seaborn as sns
import numpy as np
import pandas as pd

from Scripts.i2_parametrization.pam_parametrizer_performance_analysis import get_statistics_from_df
from Figures.Scripts.Figure2_sensitivity_error import (get_fluxomics_data,
                                                       get_reactions_to_save,
                                                       get_simulation_results_for_pams,
                                                       set_up_pams_different_parameters,
                                                       create_simulation_error_boxplot
)
from Scripts.i3_analysis.flux_kcat_distribution import create_flux_histogram_old_vs_new

N_ALT_MODELS = 8
FONTSIZE = 16

ECOLI_PHENOTYPE_DATA_PATH = os.path.join('Data', 'Ecoli_phenotypes')

MODEL_FILE_PATH = os.path.join('Models', 'iML1515.xml')

PARAM_FILE_GOTENZ = os.path.join('Results', '1_preprocessing', 'proteinAllocationModel_iML1515_EnzymaticData_250225.xlsx')
PARAM_FILE_PREPROC = os.path.join('Results', '2_parametrization',
                                     'proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx')

BEST_INDIV_RESULT_FILES = [os.path.join('Results', '2_parametrization', 'diagnostics',
                                        f'pam_parametrizer_diagnostics_{i}.xlsx') for i in range(1, N_ALT_MODELS + 1)]
PARAMETER_RESULT_FILES = [os.path.join('Results','3_analysis','parameter_files',
                                     f'proteinAllocationModel_EnzymaticData_iML1515_{i}.xlsx') for i in range(1,N_ALT_MODELS+1)]

other_colors = {'GotEnzymes': 'grey', 'After preprocessing': 'black', 'Curated':'chocolate'}
model_colors = sns.color_palette("coolwarm", n_colors=len(BEST_INDIV_RESULT_FILES))
cmap = {**{l:c for l, c in zip([f'Alternative {i + 1}' for i in range(len(BEST_INDIV_RESULT_FILES))], model_colors)},
        **other_colors}



def main_sfig1():
    #translational_sector_config
    pass

def main_sfig2():
    #compare simulation results
    pass

def main_sfig3():
    #flux distribution
    create_flux_histogram_old_vs_new([PARAM_FILE_GOTENZ,
                                      PARAM_FILE_PREPROC] + PARAMETER_RESULT_FILES,
                                     label_names=['GotEnzymes', 'After preprocessing'] \
                                                 + [f'Alternative {i}' for i in range(1, len(BEST_INDIV_RESULT_FILES) + 1)],
                                     cumulative=True,
                                     result_fig_file = os.path.join('Figures', 'SuppFig3_flux_histogram.png'),
                                     fontsize = FONTSIZE)

def main_sfig4():
    #compare error progression
    best_individual_df = None

    for i, file in enumerate(BEST_INDIV_RESULT_FILES):
        pam_param_results = pd.read_excel(file, sheet_name='Best_Individuals')
        pam_param_error = pd.read_excel(file, sheet_name='Final_Errors')
        pam_param_results = pam_param_results.merge(pam_param_error, on='run_id', how='left')
        pam_param_results['alternative'] = i + 1

        if best_individual_df is None:
            best_individual_df = pam_param_results
        else:
            best_individual_df = pd.concat([best_individual_df, pam_param_results])
    error_per_runid_config = get_statistics_from_df(best_individual_df,
                                                    group_by=['run_id'],
                                                    columns=['ga_error', 'r_squared'])
    fig = plt.figure(figsize=(16, 8))

    # R2 value progression plots
    gs = gridspec.GridSpec(1, 2, wspace=0.3)
    axs = [fig.add_subplot(gs[i]) for i in [0,1]]
    for ax, type, annotation in zip(axs,
                                    ['ga_error', 'r_squared'],
                                    ['genetic algorithm', 'PAMparametrizer']
                                    ):
        iteration = best_individual_df.drop(
            [col for col in best_individual_df.columns if col not in ['alternative', 'run_id', type]],
            axis=1).drop_duplicates(['alternative', 'run_id'], keep='last').groupby('alternative')
        for name, group in iteration:
            ax.plot(group['run_id'], group[type],
                    label=f"Alternative {name}",
                    linestyle='dashed',
                    color = cmap[f"Alternative {name}"],
                    alpha=0.5)
        ax.scatter(error_per_runid_config['run_id'], error_per_runid_config[f'{type}_mean'], color='black', label='mean')
        ax.errorbar(error_per_runid_config['run_id'], error_per_runid_config[f'{type}_mean'],
                    yerr=error_per_runid_config[f'{type}_std'], color='black', label='Mean')
        ax.set_xlabel('run id within alternative', fontsize = FONTSIZE)
        ax.set_ylabel(rf'Mean $R^{2}$ from {annotation}', fontsize = FONTSIZE)
        ax.tick_params(axis='both', which='major', labelsize=FONTSIZE)
        ax.tick_params(axis='both', which='minor', labelsize=FONTSIZE)

    axs[0].legend()
    # Add alphabet annotations
    annotations = ["A", "B"]  # Adjust based on the number of subplots
    fontsize = FONTSIZE*1.5  # Adjust as needed

    for ax, label in zip(fig.axes, annotations):
        ax.annotate(label, xy=(0, 1), xycoords="axes fraction",
                    fontsize=fontsize, fontweight='bold',
                    xytext=(-5, 5), textcoords="offset points",
                    ha="right", va="bottom")

    fig.tight_layout()
    fig.savefig(os.path.join('Figures', 'SuppFig3_pamparametrizer_performance.png'))

def main_sfig4():
    #MFA glucose comparison
    substrate_uptake_ids = ['EX_glc__D_e']

    ecoli_pams = set_up_pams_different_parameters(PARAM_FILE_GOTENZ,
                                                  PARAM_FILE_PREPROC,
                                                  BEST_INDIV_RESULT_FILES,
                                                  MODEL_FILE_PATH)
    flux_data = get_fluxomics_data(ECOLI_PHENOTYPE_DATA_PATH)
    rxns_to_save, valid_df = get_reactions_to_save(flux_data)
    substrate_rates = [
        [-f for f in flux_data.loc[substr_id, :].values] for substr_id in substrate_uptake_ids
    ]
    simulated_fluxes = get_simulation_results_for_pams(ecoli_pams,
                                                       substrate_rates=substrate_rates,
                                                       substrate_uptake_ids=substrate_uptake_ids,
                                                       rxns_to_save=rxns_to_save)
    # visualize per flux
    substrate_ids_cur = simulated_fluxes['Curated'].substrate_id
    substrate_ids = simulated_fluxes['Alternative 1'].substrate_id

    # visualize per flux
    fluxes_to_plot = rxns_to_save[3:]
    num_fluxes_on_row = int(np.ceil(len(fluxes_to_plot) / 2))
    fig, axs = plt.subplots(ncols=num_fluxes_on_row, nrows=2, figsize=[20, 8], layout='constrained')
    axs = axs.flatten()

    # upper_fig_reactions = [rxn for rxn in fluxes_to_save[:5] if rxn != 'EDA'] #EDA is not in core model
    # lower_fig_reactions = [rxn for rxn in fluxes_to_save[5:] if rxn != 'EDA']
    upper_fig_reactions = [rxn for rxn in fluxes_to_plot[:num_fluxes_on_row]]  # EDA is not in core model
    lower_fig_reactions = [rxn for rxn in rxns_to_save[num_fluxes_on_row:]]


    for alt, fluxes in simulated_fluxes.items():
        simulated_fluxes[alt] = fluxes.sort_values('substrate')
    substrate_rates = [abs(f) for f in simulated_fluxes['Curated']['substrate']]

    for i, rxn in enumerate(fluxes_to_plot):
        axs[i].scatter(valid_df['EX_glc__D_e'],valid_df[rxn])
        for alt, fluxes in simulated_fluxes.items():
            axs[i].plot(substrate_rates, [abs(f) for f in fluxes[rxn]], label=alt,
                        color=cmap[alt])
        print(flux_data)
        axs[i].set_title(flux_data.Description.loc[rxn], fontsize=FONTSIZE)
        #     axs[i].set_yticks(fontsize = fontsize)
        #     axs[i].set_xticks(fontsize = fontsize)
        axs[i].xaxis.set_tick_params(labelsize=FONTSIZE)
        axs[i].yaxis.set_tick_params(labelsize=FONTSIZE)

    handles, labels = axs[i].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.15), ncols=6, fontsize=FONTSIZE)

    fig.supylabel('Flux rate [mmol/gCDW/h]', fontsize=FONTSIZE * 1.5)
    fig.supxlabel('Glucose uptake rate [mmol/gCDW/h]', fontsize=FONTSIZE * 1.5)

    # plt.legend(loc = 'lower center')
    # plt.tight_layout()
    plt.savefig(os.path.join('Figures', 'SuppFig4_flux_comparison_glucose_intracellular.png'))

def main_sfig5():
    #MFA carbon sources
    # load exchange rates for different carbon sources by Gerosa et al. (2015) in Ecoli BW25113
    flux_csources = pd.read_excel(os.path.join(ECOLI_PHENOTYPE_DATA_PATH, 'Ecoli_phenotypes_py.xls'),
                                  sheet_name='Fluxes_Csources',
                                  engine='openpyxl',
                                  index_col=1)
    flux_csources_df = flux_csources.drop(['Flux (publication)', 'Reversibility'], axis=1)
    new_indices = []
    for i, row in flux_csources_df.iterrows():
        if isinstance(row.name, str):
            if row.name[-2:] == '_b':
                new_indices.append(row.name[:-2])
                flux_csources_df.loc[row.name] = -row
            else:
                new_indices.append(row.name)
        else:
            new_indices.append(row.name)

    flux_csources_df.index = new_indices
    flux_csources_df = flux_csources_df[flux_csources_df.index.notnull()]
    flux_csources_df = flux_csources_df.rename(
        index={'BIOMASS_Ec_iML1515_WT_75p37M': 'BIOMASS_Ecoli_core_w_GAM'}
    ).drop('Glucose (flux ratio Glc)', axis=1)
    # extract the validation data and substrate information for each carbon source
    flux_mapper = {col: flux_csources_df .index[i]
                   for i, col in enumerate(flux_csources_df .columns)}
    fluxes_to_save = []
    # Get the fluxes we want to save
    for flux in flux_csources_df .index:
        if isinstance(flux, str):
            fluxes_to_save += [f for f in flux.split(', ')]

    # parse the fluxes such that we can run and validate simulations easily
    validation_df = pd.DataFrame(columns=list(flux_csources_df .index))
    substrate_ids = []
    substrate_uptake = []
    for substrate, fluxes in flux_csources_df .items():
        substrate_ids += [flux_mapper[substrate]]
        substrate_uptake += [fluxes.loc[flux_mapper[substrate]]]
        validation_df = pd.concat([validation_df, fluxes.to_frame().T], ignore_index=True)

    validation_df.index = list(flux_mapper.values())


    ecoli_pams = set_up_pams_different_parameters(PARAM_FILE_GOTENZ,
                                                  PARAM_FILE_PREPROC,
                                                  BEST_INDIV_RESULT_FILES,
                                                  MODEL_FILE_PATH)

    simulated_fluxes = get_simulation_results_for_pams(ecoli_pams,
                                    substrate_rates=[[rate] for rate in substrate_uptake],
                                    substrate_uptake_ids=list(validation_df.index),
                                    rxns_to_save=fluxes_to_save)


    # Prepare the figure
    fig, ax = plt.subplots(figsize=(12, 6))



    create_simulation_error_boxplot(simulated_fluxes,
                                    valid_df= validation_df,
                                    fluxes_to_save=fluxes_to_save,
                                    ax=ax,fontsize = FONTSIZE,
                                    other_colors = other_colors)
    # plt.show()
    plt.savefig(os.path.join('Figures', 'SuppFig5_flux_comparison_csources.png'))

if __name__ == '__main__':
    main_sfig3()
    # main_sfig4()
    # main_sfig5()