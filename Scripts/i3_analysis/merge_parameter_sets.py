import os
from mailcap import subst

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random

from Scripts.i3_analysis.PAMparametrizer_progress_cleaned_figure import plot_valid_data, plot_simulation

from Modules.PAM_parametrizer.pam_parametrizer import PAMParametrizer
from Scripts.pam_generation_uniprot_id import set_up_ecoli_pam
from Scripts.i2_parametrization.pam_parametrizer_iML1515 import set_up_pamparametrizer

FONTSIZE = 16
FIGWIDTH = 20
FIGHEIGHT = 20
NUM_MODELS_TO_CREATE = 100
NUM_MODELS_TO_SELCT = 5
FIG_FILE_PATH = os.path.join('Results', '3_analysis', 'merged_model_performance.png')
ERROR_FILE_PATH = os.path.join('Results', '3_analysis', 'merged_models.xlsx')

def generate_randomized_enzyme_dbs(param_df_dict:dict,
                                   num_models:int = 10)-> list[pd.DataFrame]:
    # Create a dictionary to store new parameter dataframes
    randomized_dfs = []
    model = 1
    # Iterate through each dataframe in the original dict
    while model <= num_models:
        for key, df in param_df_dict.items():
            model += 1

            # Create a copy of the dataframe to modify
            new_df = df.copy()

            # Iterate over the common reactions
            for index, row in df.iterrows():
                reaction = row.rxn_id
                direction = row.direction
                enzyme = row.uniprot_id
                # Get the kcat values for this reaction across all dataframes

                kcats = [param_df_dict[df_key].loc[
                             (param_df_dict[df_key]['rxn_id'] == reaction) &
                             (param_df_dict[df_key]['direction'] == direction) &
                             ((pd.isna(enzyme)) | (df['uniprot_id'] == enzyme)), # Include enzyme condition only if it's not NaN
                             'kcat_values'
                         ].values[0] for df_key in param_df_dict]

                # Check if there are differences in kcat values
                if len(set(kcats)) > 1:
                    # Randomly select one of the differing kcat values
                    selected_kcat = random.choice(kcats)

                    # Update the kcat value in the new dataframe
                    new_df.loc[new_df['rxn_id'] == reaction, 'kcat_values'] = selected_kcat

            # Store the new dataframe in the randomized_dfs dict
            randomized_dfs.append(new_df)
    return randomized_dfs

def get_iML1515_resulting_parameter_dfs()->dict:
    # make sure all files are closed, otherwise they will not be includes
    parameter_file_dir = os.path.join('Results', '3_analysis', 'parameter_files')

    parameter_files = [os.path.join(parameter_file_dir,file) for file in os.listdir(parameter_file_dir) if 'iML1515' in file]

    param_df_dict = {}
    for file_path in parameter_files:
        try: param_df_dict[file_path] = pd.read_excel(file_path, sheet_name='ActiveEnzymes')
        except: print(file_path, ' is not an Excel file. Or is it perhaps open?')

    return param_df_dict

def initialize_parametrizer() -> set[PAMParametrizer, np.array]:
    parametrizer = set_up_pamparametrizer(-11, -0.1, kcat_increase_factor=3)
    parametrizer._init_results_objects()
    substrate_rates = sorted(parametrizer._init_validation_df([parametrizer.min_substrate_uptake_rate,
                                                               parametrizer.max_substrate_uptake_rate])['EX_glc__D_e'])
    return parametrizer, substrate_rates

def initialize_result_figure(parametrizer:PAMParametrizer,
                             substrate_rates: np.array) -> set[plt.Figure, plt.Axes]:
    fig, axs = plot_valid_data(parametrizer, fontsize=FONTSIZE, core = False)
    print('Run reference simulations')
    # fluxes = run_simulations(pamodel, substrate_rates)
    fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                     substrate_rates=substrate_rates,
                                                     sensitivity=False)
    fig, axs = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates],
                               parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
                               iteration=0, color='black')
    return fig, axs

def perform_model_simulations(parametrizer: PAMParametrizer,
                                       substrate_rates: np.array,
                                       save_fluxes: bool = True) -> list:
    fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                     substrate_rates=substrate_rates,
                                                     sensitivity=False)
    if save_fluxes:
        for flux, rate in zip(fluxes, substrate_rates):
            parametrizer.parametrization_results.add_fluxes_from_fluxdict(flux_dict=flux,
                                                                          bin_id='final',
                                                                          substrate_reaction_id= parametrizer.substrate_uptake_id,
                                                                          substrate_uptake_rate=rate,
                                                                          fluxes_abs=False)
    return fluxes

def select_best_models(error_df:pd.DataFrame, num_best_models:int = 5) -> list:
    error_df_sorted = error_df.sort_values(by='error', ascending=False)
    best_models = error_df_sorted.model.loc[:num_best_models].values
    return best_models

def plot_best_models(fig: plt.Figure, axs: plt.Axes,
                     flux_dict:dict[list], best_model_ids:list[int],
                           reactions_to_plot: list,
                           substrate_rates: np.array) -> set[plt.Figure, plt.Axes]:
    for model_id in best_model_ids:
        fluxes = flux_dict[model_id]
        fig, axs = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates],
                                   reactions_to_plot,
                                   iteration=model_id,
                                   label=f"model{model_id}", max_iteration=len(best_model_ids))
        return fig, axs


def save_simulation_plot(fig: plt.Figure, file_path: str) -> None:
    lines, labels = fig.axes[1].get_legend_handles_labels()

    fig.set_figwidth(FIGWIDTH)
    fig.set_figheight(FIGHEIGHT)
    fig.legend(lines, labels, loc='upper left', bbox_to_anchor=(0.1, 0.99), frameon=False,
               fontsize=FONTSIZE - 5)
    fig.tight_layout()
    plt.savefig(file_path)
    plt.close(fig)

if __name__ == '__main__':
    error_df = pd.DataFrame(columns=['model', 'error'])

    param_df_dict = get_iML1515_resulting_parameter_dfs()
    randomized_enzyme_dbs = generate_randomized_enzyme_dbs(param_df_dict, NUM_MODELS_TO_CREATE)

    parametrizer, substrate_rates = initialize_parametrizer()
    fig, axs = initialize_result_figure(parametrizer, substrate_rates)
    flux_dict = {}

    for i, df in enumerate(randomized_enzyme_dbs):
        #generate model with merged enzyme parameters
        pam = set_up_ecoli_pam(enzyme_db=df, sensitivity=False)
        #update the PAMparametrizer
        parametrizer.pamodel_no_sensitivity = pam

        flux_dict[i] = perform_model_simulations(parametrizer, substrate_rates, save_fluxes = True)
        error = parametrizer.calculate_final_error()
        error_df.loc[len(error_df)] = [i, error]

    best_model_nmbrs = select_best_models(error_df, NUM_MODELS_TO_SELCT)

    with pd.ExcelWriter(ERROR_FILE_PATH, mode='w') as writer:
        for best_model in best_model_nmbrs:
            randomized_enzyme_dbs[int(best_model)].to_excel(writer, sheet_name= 'model_'+best_model, index=False)

        error_df.to_excel(writer, sheet_name= 'errors')
    plot_best_models(fig, axs, flux_dict, best_model_nmbrs,
                     substrate_rates = substrate_rates,
                     reactions_to_plot = parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot)

    save_simulation_plot(fig,FIG_FILE_PATH)
    print(error_df.to_markdown())




