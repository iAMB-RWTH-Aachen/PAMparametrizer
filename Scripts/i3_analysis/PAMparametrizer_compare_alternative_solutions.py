import pandas as pd
import os
import matplotlib.pyplot as plt

from typing import Tuple

from PAModelpy.utils.pam_generation import set_up_pam
from Scripts.i2_parametrization.pam_parametrizer_iML1515 import set_up_pamparametrizer
from Scripts.i3_analysis.PAMparametrizer_progress_cleaned_figure import plot_valid_data, plot_simulation
from Modules.utils.pam_generation import create_pamodel_from_diagnostics_file



def recreate_progress_plot(best_indiv_files:list[str],
                           labels:list[str],
                           fig_file_path:str
                           ):
    ECOLI_MODEL_PATH = os.path.join('Models', 'iML1515.xml')

    FIGWIDTH = 12
    FIGHEIGHT = 12
    FONTSIZE = 20

    j=0

    parametrizer, substrate_rates = set_up_ecoli_pam_parametrizer_and_get_substrate_uptake_rates()

    substrate_rates = sorted(substrate_rates)
    fig, axs = plot_valid_data(parametrizer, fontsize=FONTSIZE)
    print('Run reference simulations')
    # fluxes = run_simulations(pamodel, substrate_rates)
    fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                                   substrate_rates = substrate_rates,
                                                                   sensitivity = False)
    fig, axs = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates]+[0],
                               parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
                               iteration=0, color='black')

    for file, label in zip(best_indiv_files, labels):
        j +=1
        print('\nAlternative ', label, ' from file ', file)
        parametrizer.pamodel = create_pamodel_from_diagnostics_file(file,
                                                                    parametrizer._pamodel.copy(copy_with_pickle = True))
        fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                                       substrate_rates=substrate_rates,
                                                                       sensitivity=False)
        fig, axs = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates]+[0],
                                   parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
                                   iteration=j + 1, max_iteration=len(best_indiv_files), label = label)


    lines, labels = fig.axes[1].get_legend_handles_labels()

    fig.set_figwidth(FIGWIDTH)
    fig.set_figheight(FIGHEIGHT)
    fig.legend(lines, labels, loc='upper left', bbox_to_anchor=(0.1, 0.99), frameon=False,
               fontsize=FONTSIZE - 5)
    fig.tight_layout()
    plt.savefig(fig_file_path)
    plt.close(fig)

def visualize_progress_plot(alternative_param_files:list[str],
                           labels:list[str],
                           fig_file_path:str
                           ):
    FIGWIDTH = 12
    FIGHEIGHT = 12
    FONTSIZE = 20

    j=0

    parametrizer, substrate_rates = set_up_ecoli_pam_parametrizer_and_get_substrate_uptake_rates()
    fig, axs = plot_valid_data(parametrizer, fontsize=FONTSIZE)
    print('Run reference simulations')
    # fluxes = run_simulations(pamodel, substrate_rates)
    fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                                   substrate_rates = substrate_rates,
                                                                   sensitivity = False)
    fig, axs = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates],
                               parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
                               iteration=0, color='black')

    for file, label in zip(alternative_param_files, labels):
        j +=1
        print('\nAlternative ', label, ' from file ', file)

        ecoli_pam = set_up_pam(file, 'Models/iML1515.xml')
        parametrizer.pamodel = ecoli_pam

        fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                                       substrate_rates=substrate_rates,
                                                                       sensitivity=False)
        fig, axs = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates],
                                   parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
                                   iteration=j + 1, max_iteration=len(alternative_param_files), label = label)


    lines, labels = fig.axes[1].get_legend_handles_labels()

    fig.set_figwidth(FIGWIDTH)
    fig.set_figheight(FIGHEIGHT)
    fig.legend(lines, labels, loc='upper left', bbox_to_anchor=(0.1, 0.99), frameon=False,
               fontsize=FONTSIZE - 5)
    fig.tight_layout()
    plt.savefig(fig_file_path)
    plt.close(fig)

def set_up_ecoli_pam_parametrizer_and_get_substrate_uptake_rates() -> Tuple:
    parametrizer = set_up_pamparametrizer(-12, -0.1, kcat_increase_factor=3)
    parametrizer._init_results_objects()
    substrate_rates = parametrizer._init_validation_df([parametrizer.min_substrate_uptake_rate,
                                                        parametrizer.max_substrate_uptake_rate])['EX_glc__D_e']
    substrate_rates = sorted(substrate_rates)
    return parametrizer, substrate_rates

if __name__ == '__main__':
    PAM_DATA_FILE_PATH = os.path.join('Data', 'proteinAllocationModel_iML1515_EnzymaticData_240730_multi.xlsx')
    # PAM_KCAT_FILES = [os.path.join('Results', '3_analysis', 'parameter_files',
    #                                f'proteinAllocationModel_EnzymaticData_iML1515_{file_nmbr}.xlsx') for file_nmbr in
    #                   range(1, 11)]
    PAM_KCAT_FILES = [os.path.join('Results', '2_parametrization', 'diagnostics',
                                   f'pam_parametrizer_diagnostics_{file_nmbr}.xlsx') for file_nmbr in
                      range(1, 11)]
    labels =  [f'alternative {i}' for i in range(1,11)]
    fig_file_path = os.path.join('Results', '3_analysis', 'pam-parametrizer_alternatives_cleaned_figure.png')

    alternative_param_files = PAM_KCAT_FILES
    recreate_progress_plot(alternative_param_files,  labels, fig_file_path)
    # visualize_progress_plot([os.path.join('Results', '3_analysis', 'parameter_files',
    #                                f'proteinAllocationModel_EnzymaticData_iML1515_1.xlsx')], ['1'], 'Results/test.png')
