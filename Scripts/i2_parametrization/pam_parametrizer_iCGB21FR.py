import os
import sys

import pandas as pd
import numpy as np
from typing import Tuple
import warnings

warnings.filterwarnings("ignore")


from Modules.PAM_parametrizer import ValidationData, HyperParameters, ParametrizationResults
from Modules.PAM_parametrizer import PAMParametrizer
from Modules.utils.pam_generation import setup_cglutanicum_pam
from PAModelpy.utils.pam_generation import increase_kcats_in_parameter_file
from Modules.utils.pamparametrizer_setup import set_up_sector_config



MAX_SUBSTRATE_UPTAKE_RATE = -0.1

def set_up_validation_data(pam_info_file: str,
                            csources: list=None) -> list[ValidationData]:
    condition2uptake = {'Gluconate': 'EX_glcn_e',
                        'Glucose': 'EX_glc__D_e',
                       'Succinate': 'EX_succ_e',
                        'Fructose': 'EX_fru_e'}
    if csources is None: csources = list(condition2uptake.keys())
    model = setup_cglutanicum_pam(pam_info_file, sensitivity =False)
    model_reactions = [rxn.id for rxn in model.reactions]

    VALID_DATA_PATH = os.path.join('Data', 'Cglutanicum_phenotypes', 'cglutanicum_phenotypes.xlsx')
    filtered_condition2uptake = {}
    for csource, uptake_rxn in condition2uptake.items():
        if (uptake_rxn in model_reactions) and (csource in csources):
            filtered_condition2uptake[csource] = uptake_rxn

    # Load the data from the sheets
    exchanges = pd.read_excel(VALID_DATA_PATH, 'Aerobic').drop(['medium','doi','Notes'], axis=1).replace(0,np.nan)

    #make validation data objects
    validation_data_objects = []
    for substrate, c_uptake_id in filtered_condition2uptake.items():
        if c_uptake_id in exchanges.columns:
            #get all rows with uptake of this carbon source and drop all empty columns
            valid_data_df = exchanges[~exchanges[c_uptake_id].isnull()].dropna(axis=1, how='all')
        if c_uptake_id not in exchanges.columns:
            print('There are no exchange rates available for', substrate)
            continue
        valid_data_df[c_uptake_id + '_ub'] = valid_data_df[c_uptake_id]

        validation_data = ValidationData(valid_data_df, c_uptake_id, [valid_data_df[c_uptake_id].min()-1, -0.1])

        #validate only exchange rates and growth rate
        validation_data._reactions_to_plot = [data for data in valid_data_df.columns if (data[-3:] != "_ub" and data!=c_uptake_id)]
        validation_data._reactions_to_validate = ['Growth', 'EX_co2_e']
        inactive_exchanges = [rxn.id for rxn in model.exchanges if rxn.id not in model.medium and rxn.id not in ['EX_lac_e']]

        validation_data.inactive_exchanges = inactive_exchanges
        if c_uptake_id == 'EX_glc__D_e':
            validation_data.sector_configs = {
                'TranslationalProteinSector':{
                'slope': model.sectors.get_by_id('TranslationalProteinSector').tps_mu[0],
                'intercept': model.sectors.get_by_id('TranslationalProteinSector').tps_0[0]
            }}
        validation_data_objects.append(validation_data)

    return validation_data_objects

def set_up_hyperparameter(processes: int,
                          gene_flow_events:int,
                          filename_extension:str,
                          num_kcats_to_mutate:int = 20,
                          threshold_iteration:int = 10):
    hyperparams = HyperParameters
    hyperparams.threshold_iteration = threshold_iteration
    hyperparams.number_of_kcats_to_mutate = num_kcats_to_mutate
    hyperparams.genetic_algorithm_filename_base = 'genetic_algorithm_run_iCGB21FR'
    hyperparams.filename_extension = filename_extension
    hyperparams.genetic_algorithm_filename_base += filename_extension

    hyperparams.genetic_algorithm_hyperparams['time_limit'] = 60000
    hyperparams.genetic_algorithm_hyperparams['processes'] = processes
    hyperparams.genetic_algorithm_hyperparams['number_gene_flow_events'] = gene_flow_events
    hyperparams.genetic_algorithm_hyperparams['number_generations'] = 5
    hyperparams.genetic_algorithm_hyperparams['error_weights'] = {'Growth': 7, 'EX_co2_e':2}
    hyperparams.genetic_algorithm_filename_base = 'genetic_algorithm_run_iCGB21FR_'
    hyperparams.genetic_algorithm_hyperparams['print_progress'] = True
    return hyperparams

def set_up_pamparametrizer(max_substrate_uptake_rate:float,
                           min_substrate_uptake_rate:float =-10,
                           pam_info_file: str = os.path.join(
                                             'Results', '1_preprocessing',
                                             'proteinAllocationModel_iCGB21FR_EnzymaticData_250915.xlsx'),
                           processes: int =4,
                           gene_flow_events: int = 4,
                           filename_extension:str = 'iCGB21FR',
                           num_kcats_to_mutate: int =20,
                           threshold_iteration:int =10,
                           c_sources:list = ['Glucose'],
                           kcat_increase_factor: int = 1):
    pam_info_file_path_out = os.path.join(
            'Results', '2_parametrization', 'proteinAllocationModel_iCGB21FR_EnzymaticData_multi.xlsx')

    increase_kcats_in_parameter_file(kcat_increase_factor,
                                         pam_info_file_path_ori=pam_info_file,
                                         pam_info_file_path_out=pam_info_file_path_out)

    pam = setup_cglutanicum_pam(pam_info_file)

    pam.GLUCOSE_EXCHANGE_RXNID = 'EX_glc__D_e'

    validation_data = set_up_validation_data(pam_info_file = pam_info_file,
                                             csources=c_sources
                                             )
    #get the min substrate uptake rate from the validation data
    for vd in validation_data:
        if vd.id == 'EX_glc__D_e': min_substrate_uptake_rate = vd.valid_data.EX_glc__D_e.sort_values(ascending=False)[0]
    hyperparameters = set_up_hyperparameter(processes,
                                            gene_flow_events,
                                            filename_extension,
                                            num_kcats_to_mutate,
                                            threshold_iteration)

    sector_configs = set_up_sector_config(pam_info_file = pam_info_file_path_out,
                                         sectors_not_related_to_growth = ['UnusedEnzymeSector', 'TranslationalProteinSector'])
    return PAMParametrizer(pamodel=pam,
                     validation_data=validation_data,
                     hyperparameters=hyperparameters,
                           sector_configs = sector_configs,
                     substrate_uptake_id='EX_glc__D_e',
                     max_substrate_uptake_rate=max_substrate_uptake_rate,
                     min_substrate_uptake_rate=min_substrate_uptake_rate)

def run_parametrizations(n_iterations:int=5,
                         pam_info_file:str = os.path.join(
                                             'Results', '1_preprocessing',
                                             'proteinAllocationModel_iCGB21FR_EnzymaticData_250523.xlsx'),
                         kcat_increase_factor: int = 1
                         ) -> None:
    for i in range(1, n_iterations+1):
        print('Working on iteration number', i, 'out of ',n_iterations)
        print('------------------------------------------------------------------------------------------------')
        pam_parametrizer = set_up_pamparametrizer(MAX_SUBSTRATE_UPTAKE_RATE,
                                                  pam_info_file=pam_info_file,
                                                  filename_extension = f'iCGB21FR_{i}',
                                                  kcat_increase_factor=kcat_increase_factor,
                                                  c_sources=['Glucose', 'Succinate', 'Fructose', 'Gluconate'])
        #
        pam_parametrizer.run(remove_subruns=True, binned='False')

if __name__ == "__main__":
    pam_info_file = os.path.join('Results', '1_preprocessing',
                                     'proteinAllocationModel_iCGB21FR_EnzymaticData_250915.xlsx')
    if len(sys.argv)>1:
        pam_info_file = sys.argv[1]

    # pam_parametrizer = set_up_pamparametrizer(MIN_SUBSTRATE_UPTAKE_RATE, MAX_SUBSTRATE_UPTAKE_RATE,
    #                      c_sources = ['Glycerol', 'Glucose', 'Succinate', 'Fructose','m-Xylene','Toluene','Benzoate', 'Octanoate'])
    # #
    # pam_parametrizer.run(remove_subruns=True, binned = 'False')
    run_parametrizations(5, pam_info_file, kcat_increase_factor=8)
# for running:
# python -m Scripts.i2_parametrization.pam_parametrizer_iML1515