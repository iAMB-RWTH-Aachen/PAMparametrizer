import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import matplotlib.ticker as ticker

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
        line = ax.plot(substrate_rates, [abs(f[r]) for f in fluxes], linewidth=5,
                       zorder=5, color=color, label = label)

    fig.canvas.draw()
    fig.canvas.flush_events()
    return fig, axs

def plot_valid_data(parametrizer, fontsize:int = 12, core = False):
    RXN_NAME_MAPPER[parametrizer.pamodel.BIOMASS_REACTION] = 'Growth rate [$h^{-1}$]'

    # plot flux changes with glucose uptake
    fig, axs = plt.subplots(2,2, dpi=100)
    valid_data = parametrizer.validation_data.get_by_id(parametrizer.substrate_uptake_id)
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