import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex

from Scripts.pam_generation import setup_ecolicore_pam
from Scripts.Testing.pam_parametrizer_ecolicore import set_up_pamparametrizer, run_simulations

RXN_NAME_MAPPER = {'EX_ac_e': 'Acetate secretion [$mmol_{ac}/g_{CDW}/h$]',
                   'EX_glc__D_e': 'Glucose uptake [$mmol_{glc}/g_{CDW}/h$]',
                   'EX_co2_e': '$CO_2$ secretion [$mmol_{CO_2}/g_{CDW}/h$]',
                   'EX_o2_e': 'Oxygen uptake [$mmol_[{O_2}/g_{CDW}/h$]',
                   'BIOMASS_Ecoli_core_w_GAM': 'Growth rate [$h^{-1}$]'}

def run_simulations(pamodel, substrate_rates) -> list:
    fluxes = []
    for substrate in substrate_rates:
        pamodel.change_reaction_bounds(rxn_id='EX_glc__D_e',
                                       lower_bound=substrate, upper_bound=0)
        print('Running simulations with ', substrate, 'mmol/g_cdw/h of substrate going into the system')
        sol_pam =pamodel.optimize()
        if pamodel.solver.status == 'optimal' and pamodel.objective.value>0:
            fluxes.append(sol_pam)
    return fluxes

def plot_simulation(fig, axs, fluxes: pd.DataFrame, substrate_rates:list, reactions_to_plot:list,
                    iteration:int = 0, color: int = None, max_iteration:int = 2) -> plt.Figure:
    if color is None:
        # adjust color to visualize progress
        # get viridis color palette
        cmap = plt.get_cmap('viridis')
        color = to_hex(cmap(iteration / (max_iteration+1)))

    label = 'Iteration ' + str(iteration)
    if iteration == 0:
        label = 'Reference'

    for r, ax in zip(reactions_to_plot, axs.flatten()):
        # plot data
        line = ax.plot(substrate_rates, [abs(f[r]) for f in fluxes], linewidth=3.5,
                       zorder=5, color=color, label = label)

    fig.canvas.draw()
    fig.canvas.flush_events()
    return fig, axs

def plot_valid_data(parametrizer, fontsize:int = 12):
    # plot flux changes with glucose uptake
    fig, axs = plt.subplots(2,2, dpi=100)

    for r, ax in zip(parametrizer.validation_data._reactions_to_plot, axs.flatten()):
        # plot data
        x = [abs(glc) for glc in parametrizer.validation_data.valid_data_df[parametrizer.substrate_uptake_id +'_ub']]
        y = [abs(data) for data in parametrizer.validation_data.valid_data_df[r]]
        ax.set_ylabel(RXN_NAME_MAPPER[r], fontsize= fontsize)
        ax.scatter(x, y,
                       color='black', marker='o', s=30, linewidths=1.3,
                       facecolors=None, zorder=0,
                       label='Literature')
        ax.set_xlabel(RXN_NAME_MAPPER[parametrizer.substrate_uptake_id], fontsize= fontsize)

    plt.tick_params(labelsize = fontsize)
    plt.ion()   # set interactive mode
    fig.tight_layout()
    fig.show()

    return fig, axs


if __name__ == '__main__':
    FIGWIDTH = 12
    FIGHEIGHT = 12
    FONTSIZE = 15

    result_file = os.path.join('Results', 'pam_parametrizer_diagnostics_ecolicore_before.xlsx')
    parametrizer = set_up_pamparametrizer(-11, -0.1)
    pamodel = parametrizer.pamodel
    step = (parametrizer.max_substrate_uptake_rate-parametrizer.min_substrate_uptake_rate)/10
    substrate_rates = list(np.arange(parametrizer.min_substrate_uptake_rate, parametrizer.max_substrate_uptake_rate, step))+[-1]

    mutated_kcats = pd.read_excel(result_file, sheet_name='Best_Individuals')

    fig, axs = plot_valid_data(parametrizer, fontsize = FONTSIZE)
    print('Run reference simulations')
    fluxes = run_simulations(pamodel, substrate_rates)
    fig, axs = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates], parametrizer.validation_data._reactions_to_plot,
                        iteration = 0,color = 'black')

    groups = mutated_kcats.groupby('run_id')
    for j, group in groups:
        print('\nRun number ', j+1)
        for i, row in group.iterrows():
            kcat_dict = {row['rxn_id']:{row['direction']:row['kcat[s-1]']}}
            pamodel.change_kcat_value(enzyme_id=row['enzyme_id'], kcats = kcat_dict)
        fluxes = run_simulations(pamodel, substrate_rates)
        fig, axs = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates],
                                   parametrizer.validation_data._reactions_to_plot,
                                   iteration = j+1, max_iteration = len(groups))

    lines, labels = fig.axes[1].get_legend_handles_labels()

    fig.set_figwidth(FIGWIDTH)
    fig.set_figheight(FIGHEIGHT)
    fig.legend(lines, labels, loc='upper left', bbox_to_anchor = (0.07,0.99), frameon =False,
               fontsize = FONTSIZE)
    fig.tight_layout()
    plt.savefig(os.path.join('Scripts', 'Results', 'pam_parametrizer_progess_cleaned_before.png'))

