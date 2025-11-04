import os
import pickle
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib as mplt
import numpy as np
from typing import Iterable, List, Union, Optional, Dict

from cobra.io import read_sbml_model
from PAModelpy.utils import set_up_pam
from PAModelpy import PAModel


def get_model_fluxes(substrate_rates: Iterable,
                     model: Union['PAModel', 'Model'])-> List[Dict[str, float]]:
    fluxes = []
    for sub in substrate_rates:
        model.reactions.EX_glc__D_e.lower_bound = -sub
        sol = model.optimize()
        fluxes.append(sol.fluxes)
    return fluxes

def change_atm_maintenance(model: Union['Model', 'PAModel'],
                           gam: Optional[float]=None,
                           ngam: Optional[float] = None,
                           biomass_rxn_id: Optional[str] = 'BIOMASS_Ec_iML1515_core_75p37M'
                           ) -> Union['PAModel', 'Model']:

    #change ngam and gam
    if ngam is not None:
        if isinstance(model, PAModel):
            model.change_reaction_bounds('ATPM', ngam)
        else:
            model.reactions.ATPM.lower_bound = ngam
        print(fr'Changing non-growth associated maintenance to {ngam} mmol_ATP/gCDW/h')
    if gam is not None:
        m_to_coeff = {model.metabolites.atp_c:-gam,
                  **{model.metabolites.get_by_id(mid): gam for mid in ['adp_c', 'pi_c', 'h_c']
                     }
                  }
        model.reactions.get_by_id(biomass_rxn_id).add_metabolites(m_to_coeff)

        print(fr'Changing growth associated maintenance to {gam} mmol_ATP/gCDW')
    return model

if __name__ == '__main__':
    NUM_MODELS = 1
    ENERGY_FIG_FILE = os.path.join('Results', '3_analysis', 'effect_of_energy_requirements.png')
    models_to_check = {f'PAM alternative {i}': os.path.join('Results', '3_analysis', 'parameter_files',
                                       f'proteinAllocationModel_EnzymaticData_iML1515_{i}.xlsx') for i in
                       range(1, NUM_MODELS + 1)}
    models_to_check = {name: set_up_pam(file, sensitivity=False) for name, file in models_to_check.items()}
    models_to_check['iML1515'] = read_sbml_model(os.path.join('Models', 'iML1515.xml'))
    colors = [mplt.colormaps['Set2'](i/4) for i in range(4)]
    substrate_rates = np.arange(0,10)

    valid_data_df = pd.read_excel(os.path.join('Data', 'Ecoli_phenotypes', 'Ecoli_phenotypes_py_rev.xls'), sheet_name='Yields')
    valid_data_df = valid_data_df[valid_data_df.Strain == 'MG1655']

    fig, axs = plt.subplots(nrows = 2, ncols = 2)
    for linestyle, (model_id, model) in zip(['--', '-', ':'], models_to_check.items()):
        fluxes = get_model_fluxes(substrate_rates = substrate_rates, model = model)
        p = pickle.dumps(model)
        model_ioj_atpm = change_atm_maintenance(model = pickle.loads(p), ngam = 3.15, gam = 53.95)
        fluxes_old_ngam = get_model_fluxes(substrate_rates = substrate_rates, model = model_ioj_atpm)

        model_more_mntnc = change_atm_maintenance(model = model_ioj_atpm, ngam = 3.15, gam = 75.55)
        fluxes_more_mntnc = get_model_fluxes(substrate_rates = substrate_rates, model = model_more_mntnc)

        model_more_mntnc1 = change_atm_maintenance(model = model_ioj_atpm, ngam = 6.86*1.1, gam = 75.55*1.1)
        fluxes_more_mntnc1 = get_model_fluxes(substrate_rates = substrate_rates, model = model_more_mntnc1)

        for ax, rxn in zip(axs.flatten(), ['EX_o2_e', 'EX_co2_e', 'BIOMASS_Ec_iML1515_core_75p37M', 'EX_ac_e']):
            ax.scatter([abs(f) for f in valid_data_df['EX_glc__D_e']], [abs(f) for f in valid_data_df[rxn]], color = 'black')
            for color, fluxlist, annotation in zip(colors,[fluxes, fluxes_old_ngam, fluxes_more_mntnc,  fluxes_more_mntnc1],
                                                       [' gam: 75.55, ngam: 6.86', ' gam: 53.95, ngam: 3.15',' gam: 75.55, ngam: 3.15',f' gam: {75.55*1.1}, ngam: {6.86*1.1}']):
                ax.plot([abs(f['EX_glc__D_e']) for f in fluxlist], [abs(f[rxn]) for f in fluxlist],
                    color = color, linestyle = linestyle, label = model_id+annotation)

            ax.set_xlabel('glc uptake rate [mmol/gcdw/h]')
            ax.set_ylabel(rxn)
    handles, labels = axs.flatten()[1].get_legend_handles_labels()

    fig.legend(labels = labels, handles = handles, loc = 'lower center',bbox_to_anchor = [0.5, 0], ncols = 2)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2)

    plt.savefig(ENERGY_FIG_FILE)


