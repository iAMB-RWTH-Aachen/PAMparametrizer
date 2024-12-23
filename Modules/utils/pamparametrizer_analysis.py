import pandas as pd
import numpy as np
import cobra

from typing import Iterable

from PAModelpy import PAModel

from Modules.utils.sector_config_functions import change_translational_sector_with_config_dict


def get_results_from_simulations(pamodel: PAModel,
                                 substrate_rates: Iterable[float],
                                 substrate_id: str = 'EX_glc__D_e',
                                 fluxes_to_save: list[str] = None,
                                 proteins_to_save:list[str] = None,
                                 transl_sector_config=True) -> dict[str, pd.DataFrame]:

    _set_up_pamodel_for_simulations(pamodel, substrate_id, transl_sector_config)
    solution_information = _set_up_solution_info(fluxes_to_save, proteins_to_save)

    for substrate in substrate_rates:
        pamodel.change_reaction_bounds(rxn_id=substrate_id,
                                       lower_bound=substrate, upper_bound=0)
        print('Running simulations with ', substrate, 'mmol/g_cdw/h of substrate going into the system')
        sol_pam = pamodel.optimize()
        if pamodel.solver.status == 'optimal' and pamodel.objective.value > 0:
            if fluxes_to_save is not None:
                solution_information['fluxes'] = save_fluxes(sol_pam, pamodel, fluxes_to_save, substrate, solution_information['fluxes'])
            if proteins_to_save is not None:
                solution_information['proteins'] = save_proteins(pamodel, proteins_to_save, substrate, solution_information['proteins'])

    return solution_information

def get_results_from_simulations_fixed_mu(pamodel: PAModel,
                                          growth_rates: Iterable[float],
                                          substrate_id: str = 'EX_glc__D_e',
                                          fluxes_to_save: list[str] = None,
                                          proteins_to_save:list[str] = None,
                                          transl_sector_config=True,
                                          method_ids: list[str]= None) -> dict[str, pd.DataFrame]:
    _set_up_pamodel_for_simulations(pamodel, substrate_id, transl_sector_config)
    solution_information = _set_up_solution_info(fluxes_to_save, proteins_to_save, method_ids)

    #leave the substrate rate open until maximal physiological accurate uptake rate
    pamodel.change_reaction_bounds(substrate_id, -15, 0)

    for method, mu in zip(method_ids, growth_rates):
        # first get the max growth rate
        pamodel.objective = pamodel.BIOMASS_REACTION
        pamodel.change_reaction_bounds(rxn_id=pamodel.BIOMASS_REACTION,
                                       lower_bound=0, upper_bound=mu)
        print('Running simulations with a growth rate of', mu, '1/h')
        pamodel.optimize()
        growth_rate = pamodel.objective.value

        #do another simulation to minimize substrate uptake rate to remove any futile cycles
        pamodel.objective = substrate_id
        pamodel.change_reaction_bounds(rxn_id=pamodel.BIOMASS_REACTION,
                                       lower_bound=growth_rate, upper_bound=growth_rate)


        sol_pam = pamodel.optimize()
        substrate_rate = pamodel.reactions.get_by_id(substrate_id).flux
        print(f'The simulated growth rate is {growth_rate} 1/h with a substrate uptake rate of {substrate_rate} mmol/g_cdw/h')
        if pamodel.solver.status == 'optimal' and growth_rate > 0:
            if fluxes_to_save is not None:
                solution_information['fluxes'] = save_fluxes(sol_pam, pamodel, fluxes_to_save, substrate_rate,
                                                             solution_information['fluxes'], method)
            if proteins_to_save is not None:
                solution_information['proteins'] = save_proteins(pamodel, proteins_to_save, substrate_rate,
                                                                 solution_information['proteins'], method)

    return solution_information

def _set_up_pamodel_for_simulations(pamodel:PAModel,
                                   substrate_id: str,
                                   transl_sector_config:bool) -> None:
    if transl_sector_config:
        transl_sector_config = {'slope': pamodel.sectors.get_by_id('TranslationalProteinSector').tps_mu[0],
                                'intercept': pamodel.sectors.get_by_id('TranslationalProteinSector').tps_0[0]}

        change_translational_sector_with_config_dict(pamodel=pamodel,
                                                     transl_sector_config = transl_sector_config,
                                                     substrate_uptake_id = substrate_id)

def _set_up_solution_info(fluxes_to_save: list[str],
                          proteins_to_save:list[str],
                          method_ids:list[str] = None
                          ) -> dict[str, pd.DataFrame]:
    solution_information = {}
    if fluxes_to_save is not None:
        solution_information['fluxes'] = pd.DataFrame(columns=['substrate'] + fluxes_to_save)
    if proteins_to_save is not None:
        solution_information['proteins'] = pd.DataFrame(
            columns=['enzyme_id', 'fraction', 'growth_rate', 'substrate_uptake'])
    if method_ids is not None:
        for key, df in solution_information.items():
            df['method'] = []
            solution_information[key]=df
    return solution_information


def save_fluxes(solution: cobra.Solution,
                pamodel:PAModel,
                fluxes_to_save: list[str],
                substrate_rate: float,
                flux_df: pd.DataFrame,
                method:str = None) -> pd.DataFrame:
    solution_flux = [substrate_rate] + [solution[rxn] if rxn in pamodel.reactions else np.NaN for rxn in
                                   fluxes_to_save]
    if method is not None: solution_flux.append(method)
    flux_df.loc[len(flux_df)] = solution_flux
    flux_df.sort_values('substrate', ascending=False)
    return flux_df

def save_proteins(pamodel:PAModel,
                proteins_to_save: list[str],
                substrate_rate: float,
                protein_df: pd.DataFrame,
                method: str = None) -> pd.DataFrame:
    total_conc = 0
    for enzid in proteins_to_save:
        enzyme = pamodel.enzymes.get_by_id(enzid)
        conc = enzyme.concentration
        total_conc += conc
        new_row = [enzid, conc, pamodel.reactions.get_by_id(pamodel.BIOMASS_REACTION).flux, substrate_rate]
        if method is not None: new_row.append(method)
        protein_df.loc[len(protein_df)] = new_row
    return protein_df

