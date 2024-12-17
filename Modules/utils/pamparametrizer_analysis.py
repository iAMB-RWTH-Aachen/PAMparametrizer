import pandas as pd
import numpy as np
import cobra

from typing import Iterable

from PAModelpy import PAModel

from Modules.utils.sector_config_functions import change_translational_sector_with_config_dict

def get_results_from_simulations(pamodel: PAModel,
                                     substrate_rates: Iterable[float],
                                     fluxes_to_save: list[str] = None,
                                     proteins_to_save:list[str] = None,
                    transl_sector_config=True) -> dict[str, pd.DataFrame]:
    if transl_sector_config:
        transl_sector_config = {'slope': pamodel.sectors.get_by_id('TranslationalProteinSector').tps_mu[0],
                                'intercept': pamodel.sectors.get_by_id('TranslationalProteinSector').tps_0[0]}

        change_translational_sector_with_config_dict(pamodel=pamodel,
                                                     transl_sector_config = transl_sector_config,
                                                     substrate_uptake_id = 'EX_glc__D_e')
    solution_information = {}
    if fluxes_to_save is not None:
        solution_information['fluxes'] = pd.DataFrame(columns=['substrate'] + fluxes_to_save)
    if proteins_to_save is not None:
        solution_information['proteins'] = pd.DataFrame(columns=['enzyme_id', 'fraction', 'growth_rate', 'substrate_uptake'])

    for substrate in substrate_rates:
        pamodel.change_reaction_bounds(rxn_id='EX_glc__D_e',
                                       lower_bound=substrate, upper_bound=0)
        print('Running simulations with ', substrate, 'mmol/g_cdw/h of substrate going into the system')
        sol_pam = pamodel.optimize()
        if pamodel.solver.status == 'optimal' and pamodel.objective.value > 0:
            if fluxes_to_save is not None:
                solution_information['fluxes'] = save_fluxes(sol_pam, pamodel, fluxes_to_save, substrate, solution_information['fluxes'])
            if proteins_to_save is not None:
                solution_information['proteins'] = save_proteins(pamodel, proteins_to_save, substrate, solution_information['proteins'])

    return solution_information

def save_fluxes(solution: cobra.Solution,
                pamodel:PAModel,
                fluxes_to_save: list[str],
                substrate_rate: float,
                flux_df: pd.DataFrame) -> pd.DataFrame:
    solution_flux = [substrate_rate] + [solution[rxn] if rxn in pamodel.reactions else np.NaN for rxn in
                                   fluxes_to_save]
    flux_df.loc[len(flux_df)] = solution_flux
    flux_df.sort_values('substrate', ascending=False)
    return flux_df

def save_proteins(pamodel:PAModel,
                proteins_to_save: list[str],
                substrate_rate: float,
                protein_df: pd.DataFrame) -> pd.DataFrame:
    total_conc = 0
    for enzid in proteins_to_save:
        enzyme = pamodel.enzymes.get_by_id(enzid)
        conc = enzyme.concentration
        total_conc += conc
        protein_df.loc[len(protein_df)] = [enzid, conc, pamodel.objective.value, substrate_rate]
    return protein_df

