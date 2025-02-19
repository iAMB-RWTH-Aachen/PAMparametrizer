import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import matplotlib.ticker as ticker
from typing import Callable

# from Scripts.i2_parametrization.pam_parametrizer_ecolicore import set_up_pamparametrizer as set_up_pamparametrizer_ecolicore
from Scripts.i2_parametrization.pam_parametrizer_iML1515 import set_up_pamparametrizer as set_up_pamparametrizer_ecoli
from Scripts.i2_parametrization.pam_parametrizer_iJN1463 import set_up_pamparametrizer as set_up_pamparametrizer_putida


FONTSIZE = 16
RXN_NAME_MAPPER = {'EX_ac_e': 'Acetate secretion [$mmol_{ac}/g_{CDW}/h$]',
                   'EX_glc__D_e': 'Glucose uptake [$mmol_{glc}/g_{CDW}/h$]',
                   'EX_co2_e': '$CO_2$ secretion [$mmol_{CO_2}/g_{CDW}/h$]',
                   'EX_o2_e': 'Oxygen uptake [$mmol_[{O_2}/g_{CDW}/h$]',
                   'BIOMASS_Ecoli_core_w_GAM': 'Growth rate [$h^{-1}$]',
                   'BIOMASS_Ec_iML1515_core_75p37M': 'Growth rate [$h^{-1}$]',
                   'BIOMASS_KT2440_WT3': 'Growth rate [$h^{-1}$]',
                   }


def plot_simulation(fig, axs,
                    fluxes: pd.DataFrame,
                    substrate_rates:list,
                    reactions_to_plot:list,
                    iteration:int = 0,
                    color: int = None,
                    max_iteration:int = 2,
                    label:str=None) -> plt.Figure:
    if color is None:
        # adjust color to visualize progress
        # get viridis color palette
        cmap = plt.get_cmap('coolwarm')
        color = to_hex(cmap(iteration / (max_iteration+1)))

    if label is None:
        label = 'Iteration ' + str(iteration)
    if iteration == 0:
        label = 'After preprocessing'

    for r, ax in zip(reactions_to_plot, axs.flatten()):
        # plot data
        line = ax.plot(substrate_rates[:-1], [abs(f[r]) for f in fluxes], linewidth=5,
                       zorder=5, color=color, label = label)

    fig.canvas.draw()
    fig.canvas.flush_events()
    return fig, axs

def plot_valid_data(parametrizer, fontsize:int = 12, core = False):
    RXN_NAME_MAPPER[parametrizer.pamodel.BIOMASS_REACTION] = 'Growth rate [$h^{-1}$]'

    # plot flux changes with glucose uptake
    fig, axs = plt.subplots(2,2, dpi=100)
    valid_data = parametrizer.validation_data.get_by_id(parametrizer.substrate_uptake_id)
    print(valid_data.valid_data)
    for r, ax in zip(valid_data._reactions_to_validate, axs.flatten()):
        # plot data
        x = [abs(glc) for glc in valid_data.valid_data[parametrizer.substrate_uptake_id + '_ub']]
        y = [abs(data) for data in valid_data.valid_data[r]]
        ax.set_ylabel(RXN_NAME_MAPPER[r], fontsize= fontsize)
        ax.scatter(x, y,
                       color='black', marker='o', s=200,
                       facecolors=None, zorder=0,
                       label='Literature')
        ax.set_xlabel(RXN_NAME_MAPPER[parametrizer.substrate_uptake_id], fontsize= fontsize)
        ax.tick_params(axis='both',labelsize = fontsize)
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))

    plt.ion()   # set interactive mode
    fig.tight_layout()
    fig.show()

    return fig, axs

def recreate_progress_plot(best_individual_df:pd.DataFrame,
                           fig_file_path:str,
                           return_error_df = False,
                           set_up_parametrizer: Callable = None):
    FIGWIDTH = 12
    FIGHEIGHT = 12
    FONTSIZE = 20

    if set_up_parametrizer is None:
        set_up_parametrizer = set_up_pamparametrizer_ecoli

    error_df = pd.DataFrame(columns=['run_id', 'error'])
    parametrizer = set_up_parametrizer(-11, -0.1, kcat_increase_factor=3)
    parametrizer._init_results_objects()
    substrate_rates = parametrizer._init_validation_df([parametrizer.min_substrate_uptake_rate,
                                                        parametrizer.max_substrate_uptake_rate])['EX_glc__D_e']
    substrate_rates = sorted(substrate_rates)
    # step = (parametrizer.max_substrate_uptake_rate - parametrizer.min_substrate_uptake_rate) / 10
    # substrate_rates = list(
    #     np.arange(parametrizer.min_substrate_uptake_rate, parametrizer.max_substrate_uptake_rate, step)) + [-1]

    fig, axs = plot_valid_data(parametrizer, fontsize=FONTSIZE)
    print('Run reference simulations')
    # fluxes = run_simulations(pamodel, substrate_rates)
    fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                                   substrate_rates = substrate_rates,
                                                                   sensitivity = False)
    fig, axs = plot_simulation(fig, axs, fluxes,
                               [abs(rate) for rate in substrate_rates]+[0],
                               parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
                               iteration=0, color='black')

    groups = best_individual_df.groupby('run_id')
    for j, group in groups:
        print('\nRun number ', j + 1)
        for i, row in group.iterrows():
            kcat_dict = {row['rxn_id']: {row['direction']: row['kcat[s-1]']}}
            parametrizer.pamodel_no_sensitivity.change_kcat_value(enzyme_id=row['enzyme_id'], kcats=kcat_dict)
        # fluxes = run_simulations(pamodel, substrate_rates)
        fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                                       substrate_rates=substrate_rates,
                                                                       sensitivity=False)
        fig, axs = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates]+[0],
                                   parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
                                   iteration=j + 1, max_iteration=len(groups))
        for flux, rate in zip(fluxes, substrate_rates):
            parametrizer.parametrization_results.add_fluxes_from_fluxdict(flux_dict=flux,
                                                                          bin_id='final',
                                                                          substrate_reaction_id= parametrizer.substrate_uptake_id,
                                                                          substrate_uptake_rate=rate,
                                                                          fluxes_abs=False)
        error = parametrizer.calculate_final_error()
        error_df.loc[len(error_df)] = [j, error]


    lines, labels = fig.axes[1].get_legend_handles_labels()

    fig.set_figwidth(FIGWIDTH)
    fig.set_figheight(FIGHEIGHT)
    fig.legend(lines, labels, loc='upper left', bbox_to_anchor=(0.1, 0.99), frameon=False,
               fontsize=FONTSIZE - 5)
    fig.tight_layout()
    plt.savefig(fig_file_path)
    plt.close(fig)
    if return_error_df: return error_df

def create_empty_plot():
    parametrizer =  set_up_pamparametrizer_ecoli(-12,-0.1)
    parametrizer._init_results_objects()
    substrate_rates = parametrizer._init_validation_df([parametrizer.min_substrate_uptake_rate,
                                                        parametrizer.max_substrate_uptake_rate])['EX_glc__D_e']
    substrate_rates = sorted(substrate_rates)
    fig, axs = plot_valid_data(parametrizer, fontsize=FONTSIZE, core =False)
    fig.set_size_inches(18.5, 10.5)
    # print('Run reference simulations')
    # # fluxes = run_simulations(pamodel, substrate_rates)
    # fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
    #                                                  substrate_rates=substrate_rates,
    #                                                  sensitivity=False)
    # fig, axs = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates],
    #                            parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
    #                            iteration=0, color='black')
    return fig, axs

def main_ecoli():
    # fig, axs = create_empty_plot()
    # plt.savefig(os.path.join('Results', 'empty_progress_plot_references.png'))
    #
    result_file = os.path.join('Results', '2_parametrization', 'diagnostics', 'pam_parametrizer_diagnostics_5.xlsx')

    best_indiv_df = pd.read_excel(result_file, sheet_name='Best_Individuals')
    #
    fig_file_path = os.path.join('Results', '3_analysis', 'pam_parametrizer_progess_cleaned_iML1515_2.png')
    recreate_progress_plot(best_indiv_df, fig_file_path)
    # result_dir = 'Results/data_reduction_results/diagnostics'
    # for file in os.listdir(result_dir):
    #     file_path = os.path.join(result_dir,file)
    #     with pd.ExcelWriter(file_path, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
    #         to_improve = pd.read_excel(file_path, sheet_name='Final_Errors')
    #         to_improve.drop_duplicates('run_id',keep='last', inplace=True)
    #
    #         # Find the last value in the 'run_id' column
    #         last_run_id = to_improve['run_id'].iloc[-1]
    #
    #         # Filter rows where 'run_id' is less than or equal to the last value
    #         filtered_df = to_improve[to_improve['run_id'] <= last_run_id]
    #
    #         filtered_df.to_excel(writer, sheet_name = 'Final_Errors', index=False)

def main_putida():
    # fig, axs = create_empty_plot()
    # plt.savefig(os.path.join('Results', 'empty_progress_plot_references.png'))
    #
    result_file = os.path.join('Results', '2_parametrization', 'diagnostics', 'pam_parametrizer_diagnostics_iJN1463_1.xlsx')

    best_indiv_df = pd.read_excel(result_file, sheet_name='Best_Individuals')
    #
    fig_file_path = os.path.join('Results', '3_analysis', 'pam_parametrizer_progess_cleaned_iJN1463_1.png')
    recreate_progress_plot(best_indiv_df,
                           fig_file_path,
                           set_up_parametrizer=set_up_pamparametrizer_putida)
    # result_dir = 'Results/data_reduction_results/diagnostics'
    # for file in os.listdir(result_dir):
    #     file_path = os.path.join(result_dir,file)
    #     with pd.ExcelWriter(file_path, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
    #         to_improve = pd.read_excel(file_path, sheet_name='Final_Errors')
    #         to_improve.drop_duplicates('run_id',keep='last', inplace=True)
    #
    #         # Find the last value in the 'run_id' column
    #         last_run_id = to_improve['run_id'].iloc[-1]
    #
    #         # Filter rows where 'run_id' is less than or equal to the last value
    #         filtered_df = to_improve[to_improve['run_id'] <= last_run_id]
    #
    #         filtered_df.to_excel(writer, sheet_name = 'Final_Errors', index=False)


if __name__ == '__main__':
    main_putida()

