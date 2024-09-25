import pandas as pd
import numpy as np
import os
from typing import Union, Callable
from PAModelpy import PAModel
import matplotlib.pyplot as plt

from Scripts.pam_generation_uniprot_id import set_up_ecoli_pam, setup_yeast_pam
from Scripts.Testing.pam_parametrizer_iML1515 import set_up_pamparametrizer
from Scripts.Testing.pam_parametrizer_yeast9 import set_up_pamparametrizer as set_up_pamparam_yeast


VALID_DATA_PATH = os.path.join('Data', 'Ecoli_phenotypes', 'Ecoli_phenotypes_py_rev.xls')
valid_data_df = pd.read_excel(VALID_DATA_PATH, sheet_name='Yields').sort_values('BIOMASS_Ec_iML1515_core_75p37M',
                                                                                    ascending = False)
max_mu = valid_data_df.at[0,'BIOMASS_Ec_iML1515_core_75p37M']

def change_kcat_with_factor(pam:PAModel,
                            multiplication_factor: int) -> None:
    for enzyme in pam.enzymes:
        kcats = enzyme.rxn2kcat.copy()
        for rxn, kcat_dict in kcats.items():
            # if the enzyme is part of a complex it has only zero kcats
            if all([val == 0 for val in kcat_dict.values()]):
                continue
            for dir, kcat in kcat_dict.items():
                kcats[rxn][dir] = kcat * multiplication_factor
        pam.change_kcat_value(enzyme.id, kcats)

def scan_kcat_factors(max_factor:int,
                      scan_figure_file_path: str,
                      min_factor:int=1,
                      stepsize:Union[float, int] = 1,
                      setup_pam_function: Callable = set_up_ecoli_pam,
                      setup_pamparametrizer_function: Callable = set_up_pamparametrizer,
                      pam_parametrizer_kwargs: dict = {'c_sources': ['Glucose'],
                                                       'threshold_iteration': 5},
                      substrate_reaction_id = 'EX_glc__D_e',
                      substrate_uptake_rates: np.arange = np.arange(-11,1,1)) -> None:

    pam_parametrizer = setup_pamparametrizer_function(substrate_uptake_rates[0],
                                                      substrate_uptake_rates[-1],
                                                      **pam_parametrizer_kwargs)
    fig, axs = pam_parametrizer.plot_valid_data()

    for mult_factor in np.arange(min_factor,max_factor+1,stepsize):
        pam_parametrizer.iteration = mult_factor

        pam = setup_pam_function(sensitivity=False)
        pam.change_reaction_bounds(substrate_reaction_id, -1e3, 0)
        change_kcat_with_factor(pam, mult_factor)

        pam_parametrizer.pamodel=pam
        fig = pam_parametrizer.plot_simulation(fig,axs, cbar_label='Scaling factor')

    plt.savefig(scan_figure_file_path)
    print("Done scanning the kcat multiplication factors")

def scan_kcat_factors_iML1515():
    scan_kcat_factors(5, os.path.join('Results', 'multifactor_scan_iML1515.png'))

def scan_kcat_factors_yeast9():
    scan_kcat_factors(max_factor=5,
                      scan_figure_file_path= os.path.join('Results', 'yeast9','multifactor_scan_yeast9.png'),
                      setup_pam_function=setup_yeast_pam,
                      setup_pamparametrizer_function=set_up_pamparam_yeast,
                      substrate_reaction_id= 'r_1714',
                      substrate_uptake_rates= np.arange(-15,1,1))

if __name__ == '__main__':
    # scan_kcat_factors_iML1515()
    scan_kcat_factors_yeast9()