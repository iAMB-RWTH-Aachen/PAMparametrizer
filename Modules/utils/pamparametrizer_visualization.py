from typing import Dict, Any
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
                   'Growth': 'Growth rate [$h^{-1}$]',
                   }

# RXN_NAME_MAPPER = {'EX_ac_e': 'Acetate secretion',
#                    'EX_glc__D_e': 'Glucose uptake',
#                    'EX_co2_e': 'CO2 secretion',
#                    'EX_o2_e': 'O2 uptake',
#                    'BIOMASS_Ecoli_core_w_GAM': 'Growth rate [$h^{-1}$]',
#                    'BIOMASS_Ec_iML1515_core_75p37M': 'Growth rate [$h^{-1}$]',
#                    'BIOMASS_KT2440_WT3': 'Growth rate [$h^{-1}$]',
#                    }


def plot_simulation(fig, axs,
                    fluxes: pd.DataFrame,
                    substrate_rates:list,
                    reactions_to_plot:list,
                    iteration:int = 0,
                    color: int = None,
                    max_iteration:int = 2,
                    label:str=None,
                    return_color: bool = False,
                    plotting_kwargs: Dict[str, Any] = {}
                    ) -> plt.Figure:
    if color is None:
        # adjust color to visualize progress
        # get viridis color palette
        cmap = plt.get_cmap('coolwarm')
        color = to_hex(cmap(iteration / (max_iteration+1)))

    if not isinstance(axs, list):axs = axs.flatten()

    if label is None:
        label = 'Iteration ' + str(iteration)
    if iteration == 0:
        label = 'After preprocessing'

    for r, ax in zip(reactions_to_plot, axs):
        # plot data
        line = ax.plot(substrate_rates, [abs(f[r]) for f in fluxes], linewidth=5,
                       zorder=5, color=color, label = label, **plotting_kwargs)

    fig.canvas.draw()
    fig.canvas.flush_events()
    if return_color: return fig, axs, color
    return fig, axs

def plot_flux_vs_experiment(ax,
                           parametrizer,
                            color,
                            fontsize):
    fluxes_dict = {}
    for substrate_id in parametrizer.substrate_uptake_ids:
        print(substrate_id)
        if substrate_id == parametrizer.substrate_uptake_id:
            substrate_range = None
        else:
            substrate_range = parametrizer.validation_data.get_by_id(substrate_id).valid_data[substrate_id]

        fluxes, substrate_range = parametrizer.run_simulations_to_plot(substrate_uptake_id=substrate_id,
                                                               substrate_rates=substrate_range,
                                                               save_fluxes_esc=False,
                                                               sensitivity=False)
        if len(fluxes) > 0:  # only plot feasible model results
            fluxes_dict[substrate_id] = fluxes

    # Plot a flux comparison to experimental data for the other fluxes
    for substr_id, fluxes in fluxes_dict.items():
        valid_data = parametrizer.validation_data.get_by_id(substr_id)
        if substr_id == parametrizer.substrate_uptake_id: continue

        feas_sampled_data = valid_data.valid_data

        for reaction in valid_data.valid_data.columns:
            if '_ub' in reaction: continue # the reaction bound is given with _ub as suffix
            exp_measurements = feas_sampled_data[reaction]
            simulations = [f[reaction] for f in fluxes]
            ax.scatter(exp_measurements, simulations, color=color)

        ax.tick_params(axis='both', labelsize=fontsize)
        ax.set_ylabel(r'Experimental flux', fontsize = fontsize)
        ax.set_xlabel(r'Simulated flux', fontsize = fontsize)

        ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useOffset=False, useMathText=False))
        ax.yaxis.get_major_formatter().set_scientific(False)  # Disable scientific notation
        ax.yaxis.get_major_formatter().set_useLocale(True)  # Ensure proper decimal formatting

        # Make the axes square
        ax.set_aspect('equal', adjustable='box')
        xlims = ax.get_xlim()
        ylims = ax.get_ylim()

        # Define the range for the diagonal
        lims = [min(xlims[0], ylims[0]), max(xlims[1], ylims[1])]
        ax.plot(lims, lims, color='black', linestyle='--', linewidth=1)
        ax.set_xlim(lims)
        ax.set_ylim(lims)

def plot_valid_data(parametrizer, axs=None, fig =None, fontsize:int = 12):
    RXN_NAME_MAPPER[parametrizer.pamodel.BIOMASS_REACTION] = 'Growth rate [$h^{-1}$]'

    # plot flux changes with glucose uptake
    if axs is None:
        fig, axs = plt.subplots(2,2, dpi=100)
        axs = axs.flatten()

    valid_data = parametrizer.validation_data.get_by_id(parametrizer.substrate_uptake_id)
    for r, ax in zip(valid_data._reactions_to_plot, axs):
        # plot data
        x = [abs(glc) for glc in valid_data.valid_data[parametrizer.substrate_uptake_id + '_ub']]
        y = [abs(data) for data in valid_data.valid_data[r]]
        ax.set_ylabel(RXN_NAME_MAPPER[r], fontsize= fontsize)
        ax.scatter(x, y,
                       color='black', marker='o', s=80,
                       facecolors=None, zorder=0,
                       label='Literature')
        if axs is None: ax.set_xlabel(RXN_NAME_MAPPER[parametrizer.substrate_uptake_id], fontsize= fontsize)
        ax.tick_params(axis='both',labelsize = fontsize)
        ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useOffset=False, useMathText=False))
        ax.yaxis.get_major_formatter().set_scientific(False)  # Disable scientific notation
        ax.yaxis.get_major_formatter().set_useLocale(True)  # Ensure proper decimal formatting

    #
    # if axs is not None:
    #     for i, ax in enumerate(axs):
    #         # Remove xticks from the first and third axes
    #         if i== in [0, 2]:
    #             ax.set_xticks([])
    #             ax.set_xticklabels([])

    plt.ion()   # set interactive mode
    # fig.tight_layout()
    fig.show()

    return fig, axs