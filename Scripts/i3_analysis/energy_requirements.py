import os
import matplotlib.pyplot as plt
import matplotlib as mplt
import numpy as np
from typing import Iterable, List, Union

from cobra.io import read_sbml_model
from PAModelpy.utils import set_up_pam


def get_model_fluxes(substrate_rates: Iterable,
                     model: Union['PAModel', 'Model'])-> List:
    fluxes
    for sub in substrate_rates:
        model.reactions.EX_glc__D_e.lower_bound = -sub
        sol = model.optimize()
        fluxes.append(sol.fluxes)
    return fluxes


if __name__ == '__main__':
    NUM_MODELS = 1
    ENERGY_FIG_FILE = os.path.join('Results', '3_analysis', 'effect_of_energy_requirements.png')
    models_to_check = {i: os.path.join('Results', '3_analysis', 'parameter_files',
                                       f'proteinAllocationModel_EnzymaticData_iML1515_{i}.xlsx') for i in
                       range(1, NUM_MODELS + 1)}
    models_to_check = {name: set_up_pam(file, sensitivity=False) for name, file in models_to_check.items()}
    models_to_check['iML1515'] = read_sbml_model(os.path.join('Models', 'iML1515.xml'))
    colors = [mplt.colormaps['Set2'](i/len(models_to_check)) for i in range(len(models_to_check))]
    substrate_rates = np.arange(0,10)

    fig, axs = plt.subplots(nrows = 2, ncols = 2)
    for color, model in zip(colors, models_to_check):
        fluxes = get_model_fluxes(substrate_rates = substrate_rates, model = model)
        fluxes_old_ngam = []
        for ax, rxn in zip(axs.flatten(), ['EX_o2_e', 'EX_co2_e', 'BIOMASS_Ec_iML1515_core_75p37M', 'EX_ac_e']):
            for fluxlist, linestyle in zip([fluxes, fluxes_old_ngam],['--', '-']):
                ax.plot([abs(f['EX_glc__D_e']) for f in fluxes], [abs(f[rxn]) for f in fluxes],
                    color = color, linestyle = linestyle)

    plt.savefig(ENERGY_FIG_FILE)


