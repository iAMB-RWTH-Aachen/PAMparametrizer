import os
import sys
import warnings
import pandas as pd
warnings.filterwarnings("ignore")

from PAModelpy.utils.pam_generation import set_up_pam, increase_kcats_in_parameter_file
from PAModelpy import Config

from Scripts.i2_parametrization.pam_parametrizer_iML1515 import (set_up_validation_data,
                                                                 set_up_hyperparameter,
                                                                 run_simulations)
from Modules.PAM_parametrizer import PAMParametrizer


MAX_SUBSTRATE_UPTAKE_RATE = -0.1
MIN_SUBSTRATE_UPTAKE_RATE = -12

def set_up_pamparametrizer(min_substrate_uptake_rate:float, max_substrate_uptake_rate: float,
                           pam_info_file:str = 'Results/1_preprocessing/proteinAllocationModel_iML1515_EnzymaticData_250219.xlsx',
                           processes: int =4,
                           gene_flow_events: int = 4,
                           filename_extension:str = 'iML1515',
                           num_kcats_to_mutate: int =10,
                           threshold_iteration:int =10,
                           c_sources:list = ['Glucose'],
                           kcat_increase_factor: int = 1)->PAMParametrizer:
    pam_info_file_path_out = os.path.join(
        'Results', '2_parametrization', 'proteinAllocationModel_mciML1515_EnzymaticData_multi.xlsx')

    increase_kcats_in_parameter_file(kcat_increase_factor,
                                     pam_info_file_path_ori=pam_info_file,
                                     pam_info_file_path_out=pam_info_file_path_out)

    config = Config()
    config.reset()

    ecoli_pam = set_up_pam(pam_info_file=pam_info_file_path_out,
                           model=os.path.join('Models', 'iML1515.xml'),
                           membrane_sector=True,
                           max_membrane_area = 0.0432987491)

    validation_data = set_up_validation_data(c_sources,
                                             pam_info_file = pam_info_file_path_out)
    hyperparameters = set_up_hyperparameter(processes,
                                            gene_flow_events,
                                            filename_extension,
                                            num_kcats_to_mutate,
                                            threshold_iteration)

    return PAMParametrizer(pamodel=ecoli_pam,
                           validation_data=validation_data,
                           hyperparameters=hyperparameters,
                           substrate_uptake_id=config.GLUCOSE_EXCHANGE_RXNID,
                           max_substrate_uptake_rate=max_substrate_uptake_rate,
                           min_substrate_uptake_rate=min_substrate_uptake_rate)

def run_multiple_parametrizations(n_iterations:int=10,
                                  pam_info_file: str = os.path.join('Scripts',
                                                                    'i2_parametrization',
                                                                    'pam_parametrizer_mcpam_iML1515.py')):
    for iteration in range(1,n_iterations+1):
        print('Working on iteration number', iteration, 'out of ', n_iterations)
        print('------------------------------------------------------------------------------------------------')
        parametrizer = set_up_pamparametrizer(MIN_SUBSTRATE_UPTAKE_RATE, MAX_SUBSTRATE_UPTAKE_RATE,
                                              pam_info_file= pam_info_file,
                                              filename_extension='mciML1515_'+str(iteration),
                                              c_sources=['Glucose'],
                                              kcat_increase_factor=3  )
        parametrizer.run(remove_subruns=True, binned = 'False')

        # need to reset best individual and computational performance df
        parametrizer.parametrization_results.best_individuals = pd.DataFrame(
            columns=['run_id', 'enzyme_id', 'direction', 'rxn_id', 'kcat[s-1]', 'ga_error'])
        parametrizer.parametrization_results.computational_time = pd.DataFrame(columns=['run_id', 'time_s', 'time_h'])

        # reset final errors for correct saving
        parametrizer.parametrization_results.final_errors = pd.DataFrame(columns=['run_id', 'r_squared'])


if __name__ == "__main__":
    pam_info_file = os.path.join(
        'Results', '1_preprocessing', 'proteinAllocationModel_iML1515_EnzymaticData_250219.xlsx')

    if len(sys.argv)>1:
        pam_info_file = sys.argv[1]

    # pam_parametrizer = set_up_pamparametrizer(MIN_SUBSTRATE_UPTAKE_RATE, MAX_SUBSTRATE_UPTAKE_RATE,
    #                                           kcat_increase_factor=3,c_sources = ['Glucose'])
    #
    # pam_parametrizer.run(remove_subruns=True, binned = 'False')
    run_multiple_parametrizations(pam_info_file=pam_info_file)
# for running:
# python -m Scripts.i2_parametrization.pam_parametrizer_iML1515