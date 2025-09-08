import os
import pandas as pd
import numpy as np
from typing import Tuple
import warnings
warnings.filterwarnings("ignore")



from PAModelpy.configuration import Config

from Modules.PAM_parametrizer import ValidationData, HyperParameters, ParametrizationResults
from Modules.PAM_parametrizer import PAMParametrizer
from Scripts.pam_generation import setup_ecolicore_pam

from Modules.utils.pamparametrizer_setup import set_up_sector_config
from PAModelpy.utils.pam_generation import increase_kcats_in_parameter_file

from Scripts.pam_generation_uniprot_id import setup_ecolicore_pam as setup_ecolicore_pam_uniprot
from Scripts.i2_parametrization.pam_parametrizer_iML1515 import set_up_hyperparameter, set_up_validation_data


# import sys
# sys.stdout = open('output.txt','wt')

MAX_SUBSTRATE_UPTAKE_RATE = -0.1
MIN_SUBSTRATE_UPTAKE_RATE = -10
config = Config()
config.reset()
RXNS_TO_VALIDATE = [config.ACETATE_EXCRETION_RXNID, config.CO2_EXHANGE_RXNID, config.OXYGEN_UPTAKE_RXNID, 'BIOMASS_Ecoli_core_w_GAM']


def set_up_validation_data_core(csources:list, pam_info_file: str) -> list[ValidationData]:
    objects = set_up_validation_data(csources=csources, pam_info_file = pam_info_file)
    for i in range(len(objects)):
        df = objects[i].valid_data.rename({'BIOMASS_Ec_iML1515_core_75p37M':'BIOMASS_Ecoli_core_w_GAM'}, axis =1)
        objects[i].valid_data = df
        objects[i].sector_configs = {}

        objects[i]._reactions_to_plot = ['BIOMASS_Ecoli_core_w_GAM', 'EX_ac_e', 'EX_o2_e', 'EX_co2_e']
    return objects

def set_up_hyperparameter() -> HyperParameters:
    hyperparams = HyperParameters
    hyperparams.threshold_iteration = 30
    hyperparams.number_of_kcats_to_mutate = 5
    hyperparams.filename_extension = 'ecolicore_false_multiple_csources2'
    hyperparams.genetic_algorithm_hyperparams['number_generations'] = 10
    hyperparams.genetic_algorithm_hyperparams['number_gene_flow_events'] = 5
    hyperparams.genetic_algorithm_filename_base = 'genetic_algorithm_run_ecolicore_'
    hyperparams.genetic_algorithm_hyperparams['print_progress'] = True
    return hyperparams

def run_simulations(pamodel, substrate_rates, rxn_to_validate = RXNS_TO_VALIDATE):
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

def set_up_pamparametrizer(min_substrate_uptake_rate:float,
                           max_substrate_uptake_rate: float,
                           pam_info_file = os.path.join('Results',
                                      '1_preprocessing',
                                      'proteinAllocationModel_iML1515_EnzymaticData_250827.xlsx'),
                           kcat_increase_factor: int = 10,
                           other_csources = False) -> PAMParametrizer:
    condition2uptake = {'Glycerol': 'EX_gly_e', 'Glucose': 'EX_glc__D_e', 'Acetate': 'EX_ac_e', 'Pyruvate': 'EX_pyr_e', 'Gluconate': 'EX_glcn_e', 'Succinate': 'EX_succ_e', 'Galactose': 'EX_gal_e', 'Fructose': 'EX_fru_e'}
    pam_info_file_path_out = os.path.join(
        'Results','2_parametrization', 'proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx')

    increase_kcats_in_parameter_file(kcat_increase_factor,
                                     pam_info_file_path_ori= pam_info_file,
                                     pam_info_file_path_out=pam_info_file_path_out)

    ecolicore_pam = setup_ecolicore_pam_uniprot(pam_info_file_path_out)
    #turn all the carbon exchanges off
    for uptake_rxn in condition2uptake.values():
        try:
            ecolicore_pam.change_reaction_bounds(uptake_rxn, lower_bound=0)
        except:
            continue
    if other_csources:
        csources = ['Glycerol', 'Glucose', 'Acetate', 'Pyruvate', 'Gluconate', 'Succinate', 'Galactose', 'Fructose']
    else:
        csources = ['Glucose']
    validation_data = set_up_validation_data_core(csources,
                                             pam_info_file='Data/proteinAllocationModel_iML1515_EnzymaticData_core.xlsx')
    hyperparameters = set_up_hyperparameter()

    sector_configs = set_up_sector_config(pam_info_file = 'Data/proteinAllocationModel_iML1515_EnzymaticData_core.xlsx',
                                         sectors_not_related_to_growth = ['UnusedEnzymeSector', 'TranslationalProteinSector'])

    return PAMParametrizer(pamodel=ecolicore_pam,
                     validation_data=validation_data,
                     hyperparameters=hyperparameters,
                        sector_configs=sector_configs,
                     substrate_uptake_id=config.GLUCOSE_EXCHANGE_RXNID,
                     max_substrate_uptake_rate=max_substrate_uptake_rate,
                     min_substrate_uptake_rate=min_substrate_uptake_rate)

if __name__ == "__main__":
    pam_parametrizer = set_up_pamparametrizer(MIN_SUBSTRATE_UPTAKE_RATE, MAX_SUBSTRATE_UPTAKE_RATE, other_csources=False)

    pam_parametrizer.run(remove_subruns=True, binned = 'False')
# for running:
# python -m Scripts.i2_parametrization.pam_parametrizer_toy_model
