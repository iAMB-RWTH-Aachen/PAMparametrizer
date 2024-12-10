import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import matplotlib.ticker as ticker

from Scripts.i2_parametrization.pam_parametrizer_ecolicore import set_up_pamparametrizer as set_up_pamparametrizer_ecolicore
from Scripts.i2_parametrization.pam_parametrizer_iML1515 import set_up_pamparametrizer as set_up_pamparametrizer_ecoli


RXN_NAME_MAPPER = {'EX_ac_e': 'Acetate secretion [$mmol_{ac}/g_{CDW}/h$]',
                   'EX_glc__D_e': 'Glucose uptake [$mmol_{glc}/g_{CDW}/h$]',
                   'EX_co2_e': '$CO_2$ secretion [$mmol_{CO_2}/g_{CDW}/h$]',
                   'EX_o2_e': 'Oxygen uptake [$mmol_[{O_2}/g_{CDW}/h$]',
                   'BIOMASS_Ecoli_core_w_GAM': 'Growth rate [$h^{-1}$]'}


def plot_simulation(fig, axs, fluxes: pd.DataFrame, substrate_rates:list, reactions_to_plot:list,
                    iteration:int = 0, color: int = None, max_iteration:int = 2, label:str=None) -> plt.Figure:
    if color is None:
        # adjust color to visualize progress
        # get viridis color palette
        cmap = plt.get_cmap('viridis')
        color = to_hex(cmap(iteration / (max_iteration+1)))

    if label is None:
        label = 'Iteration ' + str(iteration)
    if iteration == 0:
        label = 'Reference'

    for r, ax in zip(reactions_to_plot, axs.flatten()):
        # plot data
        line = ax.plot(substrate_rates[:-1], [abs(f[r]) for f in fluxes], linewidth=10,
                       zorder=5, color=color, label = label)

    fig.canvas.draw()
    fig.canvas.flush_events()
    return fig, axs

def plot_valid_data(parametrizer, fontsize:int = 12, core = False):
    if not core:
        RXN_NAME_MAPPER['BIOMASS_Ec_iML1515_core_75p37M'] = 'Growth rate [$h^{-1}$]'

    # plot flux changes with glucose uptake
    fig, axs = plt.subplots(2,2, dpi=100)
    valid_data = parametrizer.validation_data.get_by_id(parametrizer.substrate_uptake_id)
    for r, ax in zip(valid_data._reactions_to_plot, axs.flatten()):
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

def recreate_progress_plot(best_individual_df, fig_file_path, core = True, return_error_df = False):
    FIGWIDTH = 12
    FIGHEIGHT = 12
    FONTSIZE = 20

    if core:
        set_up_pamparametrizer = set_up_pamparametrizer_ecolicore
    else:
        set_up_pamparametrizer = set_up_pamparametrizer_ecoli

    error_df = pd.DataFrame(columns=['run_id', 'error'])
    parametrizer = set_up_pamparametrizer(-11, -0.1, kcat_increase_factor=3)
    parametrizer._init_results_objects()
    substrate_rates = parametrizer._init_validation_df([parametrizer.min_substrate_uptake_rate,
                                                        parametrizer.max_substrate_uptake_rate])['EX_glc__D_e']
    substrate_rates = sorted(substrate_rates)
    # step = (parametrizer.max_substrate_uptake_rate - parametrizer.min_substrate_uptake_rate) / 10
    # substrate_rates = list(
    #     np.arange(parametrizer.min_substrate_uptake_rate, parametrizer.max_substrate_uptake_rate, step)) + [-1]

    fig, axs = plot_valid_data(parametrizer, fontsize=FONTSIZE, core = core)
    print('Run reference simulations')
    # fluxes = run_simulations(pamodel, substrate_rates)
    fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                                   substrate_rates = substrate_rates,
                                                                   sensitivity = False)
    fig, axs = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates],
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
        fig, axs = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates],
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

if __name__ == '__main__':


    result_file = os.path.join('Results', 'i2_parametrization', 'diagnostics','pam_parametrizer_diagnostics_5.xlsx')
    best_indiv_df = pd.read_excel(result_file, sheet_name='Best_Individuals')

    fig_file_path = os.path.join('Scripts', 'Results', 'i3_analysis','pam_parametrizer_progess_cleaned_iML1515_2.png')
    recreate_progress_plot(best_indiv_df, fig_file_path, core = False)

