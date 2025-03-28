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
from Scripts.i3_analysis.metabolic_flux_distribution_vs_exp import main_iML1515 as plot_intracell_flux_distribution_ecoli

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
                                     result_fig_file = os.path.join('Figures', 'SuppFig4_flux_histogram.png'),
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

def main_sfig7():
    plot_intracell_flux_distribution_ecoli(fig_out=os.path.join('Figures', 'SuppFig_intracell_flux_distribution_heatmap_ecoli.png'))

if __name__ == '__main__':
    main_sfig7()
    # main_sfig4()
    # main_sfig5()