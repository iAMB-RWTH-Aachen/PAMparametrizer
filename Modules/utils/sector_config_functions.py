from typing import Iterable, Literal, Union, Tuple
from scipy.stats import linregress
import matplotlib.pyplot as plt
import pandas as pd

from PAModelpy import Config, PAModel


TransSectorConfig = dict[Literal['slope', 'intercept']]

def perform_linear_regression(x:Iterable[Union[float, int]],
                              y:Iterable[Union[float, int]]
                              )-> Tuple[float, float]:
    result = linregress(x=x,y=y)
    return result.slope, result.intercept

def reset_translational_sector(pamodel:PAModel,
                               slope:float, intercept:float,
                               new_id:str = None
                               ) -> PAModel:
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

def get_model_simulations_vs_sector(pamodel:PAModel,
                                    sub_uptake_rxn:str,
                                    rxn_id_to_relate_to:str,
                                    substrate_range: Iterable,
                                    intercept:float, slope: float,
                                    sector_name:str ='translational_protein',
                                    to_save:str=None) -> pd.DataFrame:
    """
        Runs model simulations and evaluates the relationship between a specified reaction
        and a sector variable.

        This function simulates the metabolic model over a range of substrate uptake values
        and calculates the sector variable using a linear relationship with a specified reaction.
        The results are returned as a DataFrame.

        Args:
            pamodel (PAModel): The metabolic model to be simulated.
            sub_uptake_rxn (str): The reaction ID for substrate uptake.
            rxn_id_to_relate_to (str): The reaction ID to which the sector variable
                should be related.
            substrate_range: The range of substrate uptake values to simulate.
            intercept (float): The intercept of the linear relationship used to compute
                the sector variable. For the default the units should be in g/gDW/h
            slope (float): The slope of the linear relationship used to compute
                the sector variable. For the default the units should be in g/gDW
            sector_name (str, optional): The name of the sector variable. Defaults to 'translational_protein'.
            to_save (str, optional): An additional column to include in the output DataFrame. Defaults to None.

        Returns:
            pd.DataFrame: A DataFrame containing the simulation results with relevant columns
            based on the specified reaction and sector variable.

        Example:
            ```python
            results = get_model_simulations_vs_sector(
                pamodel=my_model,
                sub_uptake_rxn="EX_glc__D_e",
                rxn_id_to_relate_to="BIOMASS_Ecoli",
                substrate_range=np.linspace(0, 10, 50),
                intercept=1.0,
                slope=0.5,
                sector_name="translational_protein"
            )
            ```
        """
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

def run_simulations(pamodel:PAModel,
                    substrate_rates:Iterable[float],
                    sub_uptake_id = 'EX_glc__D_e'
                    ) -> list:
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


def plot_translational_protein_vs_mu(literature:pd.DataFrame,
                                     results:pd.DataFrame,
                                     protein_fraction:float,
                                     measured_protein_fraction:float,
                                     oxygen_results:pd.DataFrame=None,
                                     oxygen_rxn_id:str=None,
                                     return_fig:bool = False,
                                     configuration:Config = None,
                                     literature_label:str = 'Schmidt et al (2016)',
                                     model_label:str = 'new iML1515 PAM'
                                     )->None:

    if configuration is None:
        configuration  = Config().reset()
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plotting the literature data
    ax.scatter(
        literature['Experimental Growth rate'],
        literature['Translation, ribosomal structure and biogenesis'] / 100 * measured_protein_fraction / protein_fraction,
        label=literature_label, color='blue', s=50, zorder=5
    )

    # Plotting the no limitation line
    ax.plot(
        results[configuration.BIOMASS_REACTION],
        abs(results['translational_protein']),
        label=model_label, linewidth=5, color='red', zorder=10
    )

    # Plotting oxygen limitation results if provided
    if oxygen_results is not None:
        for o2_level, result in oxygen_results.groupby(oxygen_rxn_id):
            ax.plot(
                result[configuration.BIOMASS_REACTION],
                abs(result['translational_protein']),
                linewidth=5,
                label=f'{round(o2_level * 100) / 100} $mmol_{{O_{{2}}}}/g_{{CDW}}/h$'
            )

    # Setting axes labels
    ax.set_xlabel('Growth rate (h$^{-1}$)', fontsize=12)
    ax.set_ylabel('Translational protein sector (g$_{translationalprotein}$/g$_{protein}$)', fontsize=12)

    # Adding a legend and setting the layout
    ax.legend()
    # ax.grid(True)
    ax.set_title('Glucose Limitation', fontsize=14)
    fig.tight_layout()

    if return_fig: return fig, ax
    plt.show()


def plot_unused_protein_vs_mu(results:pd.DataFrame,
                              biomass_rxn:str
                              )->None:
    # plotting
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plotting the no limitation line
    ax.plot(
        results[biomass_rxn],
        abs(results['unused_protein']),
        label='Simulations', linewidth=5, color='red', zorder=10
    )

    # Setting axes labels
    ax.set_xlabel('Growth rate (h$^{-1}$)', fontsize=12)
    ax.set_ylabel('Unused enzyme sector ($g_{unusedprotein}/g_{CDW}$)', fontsize=12)


    ax.set_title('Substrate Limitation', fontsize=14)
    fig.tight_layout()

    plt.show()


def change_translational_sector_with_config_dict(pamodel:PAModel,
                                                 transl_sector_config:dict,
                                                 substrate_uptake_id:str
                                                 ) -> None:
    pamodel.constraints[pamodel.TOTAL_PROTEIN_CONSTRAINT_ID].lb = 0 #need to set the lb to 0 to prevent errors in the setter methods

    pamodel.change_sector_parameters(pamodel.sectors.get_by_id('TranslationalProteinSector'),
                                              slope=transl_sector_config['slope'],
                                              intercept=transl_sector_config['intercept'],
                                    lin_rxn_id=substrate_uptake_id
                                     )

    # pamodel.constraints[pamodel.TOTAL_PROTEIN_CONSTRAINT_ID].lb = pamodel.constraints[pamodel.TOTAL_PROTEIN_CONSTRAINT_ID].ub #reset


def get_translational_sector_config(pamodel:PAModel,
                                    substrate_id: str,
                                    substrate_range:Iterable[Union[int, float]],
                                    rxn_to_relate_to: str = None
                                    )->TransSectorConfig:
    #generate a pam with only the translational sector
    #make a copy of the pam using pickle (the copy method for some reason does not work properly)
    pamtransl = pamodel.copy(copy_with_pickle=True)

    #remove total protein to remove protein relations
    pamtransl.remove_cons_vars([pamtransl.constraints[pamtransl.TOTAL_PROTEIN_CONSTRAINT_ID]])
    trans_sector = pamodel.sectors.TranslationalProteinSector

    #get the simulations for relatively low substrate uptake rates
    if rxn_to_relate_to is None: rxn_to_relate_to = pamodel.BIOMASS_REACTION
    simulation_results = get_model_simulations_vs_sector(pamodel = pamtransl,
                                                             sub_uptake_rxn = substrate_id,
                                                             rxn_id_to_relate_to = rxn_to_relate_to,
                                                             substrate_range = substrate_range,
                                                             intercept = trans_sector.tps_0[0],
                                                             slope = trans_sector.tps_mu[0])
    if len(simulation_results)==0 or len(simulation_results[substrate_id].unique()) ==1: # model was infeasible for all datapoints or did not change at all
        slope, intercept = trans_sector.tps_mu[0],trans_sector.tps_0[0],
    else:
        slope, intercept = perform_linear_regression(
            x=simulation_results[substrate_id], y=simulation_results["translational_protein"])
            # slope, intercept = slope * self.MEASURED_PROTEIN_FRACTION, intercept*self.MEASURED_PROTEIN_FRACTION
    return {"slope":slope,"intercept":intercept}