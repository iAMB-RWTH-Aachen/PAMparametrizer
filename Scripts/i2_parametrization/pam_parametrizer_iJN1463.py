import os
import pandas as pd
import numpy as np
from typing import Tuple
import warnings
warnings.filterwarnings("ignore")



from PAModelpy.configuration import Config

from Modules.PAM_parametrizer import ValidationData, HyperParameters, ParametrizationResults
from Modules.PAM_parametrizer import PAMParametrizer
from Scripts.pam_generation_uniprot_id import setup_pputida_pam
from PAModelpy.utils.pam_generation import increase_kcats_in_parameter_file



# import sys
# sys.stdout = open('output.txt','wt')

MAX_SUBSTRATE_UPTAKE_RATE = -0.1
MIN_SUBSTRATE_UPTAKE_RATE = -20
config = Config()
config.reset()

def set_up_validation_data(csources: list=None) -> list[ValidationData]:
    condition2uptake = {'Glycerol': 'EX_glyc_e', 'Glucose': 'EX_glc__D_e',
                        'Octanoate': 'EX_octa_e', 'Toluene': 'TOLtex', 'm-Xylene': 'M_Xylt1',
                         'Succinate': 'EX_succ_e', 'Benzoate': 'EX_bz_e',
                        'Fructose': 'EX_fru_e'}
    if csources is None: csources = list(condition2uptake.keys())
    model = setup_pputida_pam(sensitivity =False)
    model_reactions = [rxn.id for rxn in model.reactions]

    VALID_DATA_PATH = os.path.join('Data', 'Pputida_phenotypes', 'pputida_phenotypes.xlsx')
    filtered_condition2uptake = {}
    for csource, uptake_rxn in condition2uptake.items():
        if (uptake_rxn in model_reactions) and (csource in csources):
            filtered_condition2uptake[csource] = uptake_rxn

    # Load the data from the sheets
    exchanges = pd.read_excel(VALID_DATA_PATH, 'exchanges')
    exchanges = exchanges[exchanges.Strain == "KT2440"].drop(['Strain', 'Reference', 'Info'], axis=1).dropna(axis=1, how='all')

    # get netto glucose uptake
    exchanges['EX_glc__D_e'] = exchanges['EX_glc__D_e'] + exchanges['EX_glcn_e'] + exchanges['EX_25dkglcn_e']

    # Dictionary to map carbon sources to exchange reactions and fluxomics sheets
    # not fructose, because that is possibly mixotrophic growth
    fluxomics_csources = {}
    for rxn_id, sheet_name in zip(["EX_bz_e", "EX_glc__D_e"],["fluxomics_benzoate","fluxomics_glucose"]):
        df = pd.read_excel(VALID_DATA_PATH, sheet_name)
        df = df[df.Strain == "KT2440"].drop('Strain', axis=1)
        fluxomics_csources[rxn_id] = df

    #make validation data objects
    validation_data_objects = []
    for substrate, c_uptake_id in filtered_condition2uptake.items():
        if c_uptake_id in exchanges.columns:
            #get all rows with uptake of this carbon source and drop all empty columns
            valid_data_df = exchanges[~exchanges[c_uptake_id].isnull()].dropna(axis=1, how='all')
        if c_uptake_id in fluxomics_csources.keys():
            fluxomics_data = fluxomics_csources[c_uptake_id].drop('Reference', axis=1)
            if c_uptake_id in exchanges.columns:
                valid_data_df = pd.concat([fluxomics_data, valid_data_df])
            else:
                valid_data_df = fluxomics_data
        if (c_uptake_id not in fluxomics_csources.keys()) and (c_uptake_id not in exchanges.columns):
            print('There are no exchange rates available for', substrate)
            continue
        valid_data_df[c_uptake_id + '_ub'] = valid_data_df[c_uptake_id]

        validation_data = ValidationData(valid_data_df, c_uptake_id, [valid_data_df[c_uptake_id].min()-1, -0.1])
        #validate only exchange rates and growth rate

        validation_data._reactions_to_plot = [data for data in valid_data_df.columns if data[-3:] != "_ub"]
        validation_data._reactions_to_validate = [col for col in valid_data_df.columns if
                                                  ('EX_' in col) and (col[-3:] != "_ub")] + ['BIOMASS_KT2440_WT3']

        if c_uptake_id == 'EX_glc__D_e':
            validation_data.translational_sector_config = {
                'slope': model.sectors.get_by_id('TranslationalProteinSector').tps_mu[0],
                'intercept': model.sectors.get_by_id('TranslationalProteinSector').tps_0[0]
            }

            validation_data._reactions_to_plot = ['BIOMASS_KT2440_WT3', 'EDD','MDH', 'EX_glcn_e', 'EX_25dkglcn_e']
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
    hyperparams.genetic_algorithm_filename_base = 'genetic_algorithm_run_iML1515'
    hyperparams.filename_extension = filename_extension
    hyperparams.genetic_algorithm_filename_base += filename_extension

    hyperparams.genetic_algorithm_hyperparams['time_limit'] = 60000
    hyperparams.genetic_algorithm_hyperparams['processes'] = processes
    hyperparams.genetic_algorithm_hyperparams['number_gene_flow_events'] = gene_flow_events
    hyperparams.genetic_algorithm_hyperparams['number_generations'] = 5
    hyperparams.genetic_algorithm_filename_base = 'genetic_algorithm_run_iJN1463_'
    hyperparams.genetic_algorithm_hyperparams['print_progress'] = True
    return hyperparams

def run_simulations(pamodel, substrate_rates, rxn_to_validate):
    result_df = pd.DataFrame(columns= rxn_to_validate)

    for substrate in substrate_rates:
        pamodel.change_reaction_bounds(rxn_id=config.GLUCOSE_EXCHANGE_RXNID,
                                       lower_bound=substrate, upper_bound=0)
        print('Running simulations with ', substrate, 'mmol/g_cdw/h of substrate going into the system')
        pamodel.optimize()
        if pamodel.solver.status == 'optimal' and pamodel.objective.value>0:
            results_row = []
            for rxn in rxn_to_validate:
                results_row += [pamodel.reactions.get_by_id(rxn).flux]

            result_df.loc[len(result_df)] = [substrate] + results_row
    return result_df

def set_up_pamparametrizer(min_substrate_uptake_rate:float, max_substrate_uptake_rate: float,
                           processes: int =4,
                           gene_flow_events: int = 4,
                           filename_extension:str = 'iJN1463',
                           num_kcats_to_mutate: int =20,
                           threshold_iteration:int =10,
                           c_sources:list = ['Glucose'],
                           kcat_increase_factor: int = 1):
    if kcat_increase_factor != 1:
        pam_info_file_path_out = os.path.join(
            'Results', '2_parametrization', 'proteinAllocationModel_iJN1463_EnzymaticData_multi.xlsx')

        increase_kcats_in_parameter_file(kcat_increase_factor,
                                         pam_info_file_path_ori=os.path.join(
                                             'Results', '1_preprocessing',
                                             'proteinAllocationModel_iJN1463_EnzymaticData_250117.xlsx'),
                                         pam_info_file_path_out=pam_info_file_path_out)

    pputida_pam = setup_pputida_pam()
    #close off all exchanges not necessary for medium, as pputida doesn't exchange anything
    for exchange in pputida_pam.exchanges:
        if exchange.id not in pputida_pam.medium:
            pputida_pam.change_reaction_bounds(exchange.id, 0,0)

    pputida_pam.GLUCOSE_EXCHANGE_RXNID = 'EX_glc__D_e'


    validation_data = set_up_validation_data(c_sources)
    hyperparameters = set_up_hyperparameter(processes,
                                            gene_flow_events,
                                            filename_extension,
                                            num_kcats_to_mutate,
                                            threshold_iteration)
    return PAMParametrizer(pamodel=pputida_pam,
                     validation_data=validation_data,
                     hyperparameters=hyperparameters,
                     substrate_uptake_id='EX_glc__D_e',
                     max_substrate_uptake_rate=max_substrate_uptake_rate,
                     min_substrate_uptake_rate=min_substrate_uptake_rate)

def run_parametrizations(n_iterations:int=5) -> None:
    for i in range(1, n_iterations+1):
        print('Working on iteration number', i, 'out of ',n_iterations)
        print('------------------------------------------------------------------------------------------------')
        pam_parametrizer = set_up_pamparametrizer(MIN_SUBSTRATE_UPTAKE_RATE, MAX_SUBSTRATE_UPTAKE_RATE,
                                                  filename_extension = f'iJN1463_{i}',
                                                  c_sources=['Glycerol', 'Glucose', 'Succinate', 'Fructose', 'm-Xylene',
                                                             'Toluene', 'Benzoate', 'Octanoate'])
        #
        pam_parametrizer.run(remove_subruns=True, binned='False')

if __name__ == "__main__":
    # pam_parametrizer = set_up_pamparametrizer(MIN_SUBSTRATE_UPTAKE_RATE, MAX_SUBSTRATE_UPTAKE_RATE,
    #                      c_sources = ['Glycerol', 'Glucose', 'Succinate', 'Fructose','m-Xylene','Toluene','Benzoate', 'Octanoate'])
    # #
    # pam_parametrizer.run(remove_subruns=True, binned = 'False')
    run_parametrizations()
# for running:
# python -m Scripts.i2_parametrization.pam_parametrizer_iML1515