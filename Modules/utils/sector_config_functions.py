from scipy.stats import linregress
# import plotly.graph_objs as go
import matplotlib.pyplot as plt
import pandas as pd

def perform_linear_regression(x, y):
    result = linregress(x=x,y=y)
    return result.slope, result.intercept

def reset_translational_sector(pamodel, slope, intercept, new_id = None):
    transl_sector = pamodel.sectors.get_by_id('TranslationalProteinSector')
    pamodel.remove_sectors([transl_sector])
    if new_id is not None:
        transl_sector.id_list = [new_id]
    transl_sector.tps_mu = [slope]
    transl_sector.tps_0 = [intercept]
    transl_sector.set_slope()
    transl_sector.set_intercept()
    pamodel.add_sectors([transl_sector])
    return pamodel

def get_model_simulations_vs_sector(pamodel, sub_uptake_rxn, rxn_id_to_relate_to,
                                    substrate_range, intercept, slope,
                                    sector_name ='translational_protein',
                                    to_save:str=None):
    # run the model
    simulations_transl_vs_subst = run_simulations(pamodel, substrate_range, sub_uptake_rxn)
    if sub_uptake_rxn == rxn_id_to_relate_to:
        columns = [sub_uptake_rxn, sector_name]

    else:
        columns = [sub_uptake_rxn, rxn_id_to_relate_to, sector_name]
    if to_save is not None: columns.append(to_save)

    simulation_results = pd.DataFrame(columns=columns)
    for i, sim in enumerate(simulations_transl_vs_subst):
        sim[sector_name] = intercept + slope * sim[rxn_id_to_relate_to]
        simulation_results.loc[i] = [sim[col] for col in columns]

    return simulation_results

def run_simulations(pamodel, substrate_rates, sub_uptake_id = 'EX_glc__D_e') -> list:
    fluxes = []
    for substrate in substrate_rates:
        if substrate<0:
            pamodel.change_reaction_bounds(rxn_id=sub_uptake_id,
                                       lower_bound=substrate, upper_bound=0)
        else:
            pamodel.change_reaction_bounds(rxn_id=sub_uptake_id,
                                       lower_bound=0, upper_bound=substrate)

        print('Running simulations with ', substrate, 'mmol/g_cdw/h of substrate going into the system')
        sol_pam =pamodel.optimize()
        if pamodel.solver.status == 'optimal' and pamodel.objective.value>0:
            fluxes.append(sol_pam.fluxes)
    return fluxes


def plot_translational_protein_vs_mu(literature, results,
                                     protein_fraction, measured_protein_fraction,
                                     oxygen_results=None, oxygen_rxn_id=None, return_fig = False):
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plotting the literature data
    ax.scatter(
        literature['Experimental Growth rate'],
        literature['Translation, ribosomal structure and biogenesis'] / 100 * measured_protein_fraction / protein_fraction,
        label='Schmidt et al (2016)', color='blue', s=50, zorder=5
    )

    # Plotting the no limitation line
    ax.plot(
        results['BIOMASS_Ec_iML1515_core_75p37M'],
        abs(results['translational_protein']),
        label='new iML1515 PAM', linewidth=5, color='red', zorder=10
    )

    # Plotting oxygen limitation results if provided
    if oxygen_results is not None:
        for o2_level, result in oxygen_results.groupby(oxygen_rxn_id):
            ax.plot(
                result['BIOMASS_Ec_iML1515_core_75p37M'],
                abs(result['translational_protein']),
                linewidth=5,
                label=f'{round(o2_level * 100) / 100} $mmol_{{O_{{2}}}}/g_{{CDW}}/h$'
            )

    # Setting axes labels
    ax.set_xlabel('Growth rate (h$^{-1}$)', fontsize=12)
    ax.set_ylabel('Translational protein sector (g$_{translationalprotein}$/g$_{protein}$)', fontsize=12)

    # # Adding annotations for intercept and slope
    # ax.annotate(f"Intercept: {round(interc * 1000) / 1000} $g_{{tprot}}/mmol_{{s}}$",
    #             xy=(0.34, 0.5), xycoords='axes fraction', fontsize=10, ha='left')
    # ax.annotate(f"Slope: {round(slope * 1000) / 1000} $g_{{tprot}}/mmol_{{s}}/h$",
    #             xy=(0.30, 0.47), xycoords='axes fraction', fontsize=10, ha='left')

    # Adding a legend and setting the layout
    ax.legend()
    # ax.grid(True)
    ax.set_title('Glucose Limitation', fontsize=14)
    fig.tight_layout()

    if return_fig: return fig, ax
    plt.show()


def plot_unused_protein_vs_mu(results, ups_molmass, biomass_rxn):
    # plotting
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plotting the no limitation line
    ax.plot(
        results['BIOMASS_Ec_iML1515_core_75p37M'],
        abs(results['translational_protein']),
        label='Simulations', linewidth=5, color='red', zorder=10
    )

    # Setting axes labels
    ax.set_xlabel('Growth rate (h$^{-1}$)', fontsize=12)
    ax.set_ylabel('Unused enzyme sector ($g_{unusedprotein}/g_{CDW}$)', fontsize=12)
    ax.set_ylim([0,0.2])


    ax.set_title('Glucose Limitation', fontsize=14)
    fig.tight_layout()

    plt.show()


    # fig = go.Figure()
    # fig.add_trace(go.Scatter(x=results[biomass_rxn],
    #                          y=results['unused_protein'],
    #                          mode='lines', line=dict(width=5)))
    # fig.update_xaxes(title_text='Growth rate ($h^{-1}$)', range =[0,1])
    # fig.update_yaxes(title_text='Unused enzyme sector ($g_{unusedprotein}/g_{CDW}$)', range =[0,0.2])
    #
    # fig.update_traces(marker_size=10)
    # fig.update_layout(font=dict(size=12), title='Glucose limitation')
    # fig.show()

def change_translational_sector_with_config_dict(pamodel, transl_sector_config:dict, substrate_uptake_id:str) -> None:
    pamodel.change_sector_parameters(pamodel.sectors.get_by_id('TranslationalProteinSector'),
                                              slope=transl_sector_config['slope'],
                                              intercept=transl_sector_config['intercept'],
                                    lin_rxn_id=substrate_uptake_id
                                     )