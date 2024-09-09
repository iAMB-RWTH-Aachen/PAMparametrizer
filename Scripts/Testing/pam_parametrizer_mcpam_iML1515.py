import os
import pandas as pd
import numpy as np
from typing import Tuple
import warnings
warnings.filterwarnings("ignore")



from PAModelpy.configuration import Config

from Modules.PAM_parametrizer import ValidationData, HyperParameters, ParametrizationResults
from Modules.PAM_parametrizer import PAMParametrizer
from Scripts.pam_generation import setup_ecoli_pam
from Scripts.mcpam_generation_uniprot_id import set_up_ecoli_pam as setup_ecoli_pam_uniprot


# import sys
# sys.stdout = open('output.txt','wt')

MAX_SUBSTRATE_UPTAKE_RATE = -0.1
MIN_SUBSTRATE_UPTAKE_RATE = -10
config = Config()
config.reset()
RXNS_TO_VALIDATE = [config.ACETATE_EXCRETION_RXNID,  config.OXYGEN_UPTAKE_RXNID, config.BIOMASS_REACTION, config.CO2_EXHANGE_RXNID]


def set_up_validation_data(csources: list) -> list[ValidationData]:
    condition2uptake = {'Glycerol': 'EX_gly_e', 'Glucose': 'EX_glc__D_e', 'Acetate': 'EX_ac_e', 'Pyruvate': 'EX_pyr_e',
                        'Gluconate': 'EX_glcn_e', 'Succinate': 'EX_succ_e', 'Galactose': 'EX_gal_e',
                        'Fructose': 'EX_fru_e'}
    model = setup_ecoli_pam_uniprot(pam_info_file= os.path.join('Data', 'mcPAM_iML1515_EnzymaticData_new.xlsx'))
    model_reactions = [rxn.id for rxn in model.reactions]

    DATA_DIR = os.path.join(os.getcwd(), 'Data')
    VALID_DATA_PATH = os.path.join(DATA_DIR, 'Ecoli_phenotypes', 'Ecoli_phenotypes_py_rev.xls')
    valid_data_csources, condition2uptake = get_validation_data_df_other_csources(condition2uptake, model_reactions)

    validation_data_objects = []

    for csource in csources:
        if csource == 'Glucose':
            validation_data = set_up_valid_data_glucose(VALID_DATA_PATH)
            validation_data.translational_sector_config = {
                'slope': model.sectors.get_by_id('TranslationalProteinSector').tps_mu[0],
                'intercept': model.sectors.get_by_id('TranslationalProteinSector').tps_0[0]
            }
            validation_data_objects.append(validation_data)
        elif csource in condition2uptake.keys():
            validation_data = set_up_valid_data_csource_not_glucose(
                valid_data_csources, csource, condition2uptake)
            validation_data_objects.append(validation_data)
    return validation_data_objects


def get_validation_data_df_other_csources(condition2uptake: dict, model_reactions: list) -> Tuple[pd.DataFrame, dict]:
    VALID_DATA_CSOURCES_PATH = os.path.join('Data', 'Ecoli_phenotypes', 'fluxomics_datasets_gerosa.xlsx')
    valid_data_csources = pd.read_excel(VALID_DATA_CSOURCES_PATH, sheet_name='Gerosa et al')

    # only retain those carbon sources which are take up in the model
    filtered_condition2uptake = {}
    for csource, uptake_rxn in condition2uptake.items():
        if uptake_rxn in model_reactions:
            filtered_condition2uptake[csource] = uptake_rxn

    valid_data_csources = valid_data_csources[valid_data_csources['condition'].isin(filtered_condition2uptake.keys())]
    # only retain those reactions which are in the model
    valid_data_csources = valid_data_csources[valid_data_csources['reaction'].isin(list(model_reactions))]
    return valid_data_csources, filtered_condition2uptake

def set_up_valid_data_glucose(file_path: str) -> ValidationData:
    valid_data_df = pd.read_excel(file_path, sheet_name='Yields')
    valid_data_df = valid_data_df.rename(columns={config.GLUCOSE_EXCHANGE_RXNID: config.GLUCOSE_EXCHANGE_RXNID + '_ub'})
    # valid_data_df = valid_data_df[(valid_data_df.Reference != 'Folsom 2015') & (
    #             valid_data_df.Reference != 'Fischer 2003')]  # valid_data_df.Reference != 'Folsom 2015') &
    validation_data = ValidationData(valid_data_df, config.GLUCOSE_EXCHANGE_RXNID, [MIN_SUBSTRATE_UPTAKE_RATE, MAX_SUBSTRATE_UPTAKE_RATE])
    validation_data._reactions_to_plot = RXNS_TO_VALIDATE
    validation_data._reactions_to_validate = RXNS_TO_VALIDATE
    return validation_data

def set_up_valid_data_csource_not_glucose(valid_data_csources: pd.DataFrame, csource:str,
                                          condition2uptake: dict) -> ValidationData:
    valid_data_df = valid_data_csources[valid_data_csources.condition == csource]
    valid_data_df = valid_data_df[['reaction', 'measured']].set_index('reaction').T

    valid_data_df[condition2uptake[csource] + '_ub'] = valid_data_df[condition2uptake[csource]]
    validation_data = ValidationData(valid_data_df, condition2uptake[csource], [-30, 0])
    validation_data._reactions_to_plot = list(valid_data_df.columns)
    validation_data._reactions_to_validate = [col for col in valid_data_df.columns if not '_ub' in col]
    return validation_data

def set_up_hyperparameter(processes: int,
                          gene_flow_events:int,
                          filename_extension:str,
                          num_kcats_to_mutate:int = 4):
    hyperparams = HyperParameters
    hyperparams.threshold_iteration = 10
    hyperparams.number_of_kcats_to_mutate = num_kcats_to_mutate
    hyperparams.genetic_algorithm_filename_base = 'genetic_algorithm_run_mcecoli_'
    hyperparams.filename_extension = filename_extension
    hyperparams.genetic_algorithm_filename_base += filename_extension

    hyperparams.genetic_algorithm_hyperparams['processes'] = processes
    hyperparams.genetic_algorithm_hyperparams['number_gene_flow_events'] = gene_flow_events
    hyperparams.genetic_algorithm_hyperparams['number_generations'] = 6
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

def set_up_pamparametrizer(min_substrate_uptake_rate:float, max_substrate_uptake_rate: float,
                           processes: int =2,
                           gene_flow_events: int = 2,
                           filename_extension:str = 'iML1515',
                           num_kcats_to_mutate: int =4,
                           c_sources:list = ['Glucose']):
    pam_file = pd.read_excel(os.path.join('Data', 'mcPAM_iML1515_EnzymaticData.xlsx'), sheet_name = 'ActiveEnzymes')
    pam_file = pam_file.melt(id_vars=[col for col in pam_file.columns if col not in ['kcat_b', 'kcat_f']],
                    value_vars=['kcat_b', 'kcat_f'],
                    var_name='direction',
                    value_name='kcat_values')
    pam_file['direction'] = pam_file['direction'].replace({'kcat_b': 'b', 'kcat_f': 'f'})
    pam_file['kcat_values'] = pam_file['kcat_values'].fillna(11)
    pam_file['kcat_values'].loc[pam_file['kcat_values'] == 175] = 2000
    with pd.ExcelWriter(os.path.join('Data', 'mcPAM_iML1515_EnzymaticData_new.xlsx'), mode='a', engine='openpyxl',if_sheet_exists='replace') as writer:
        pam_file.to_excel(writer, sheet_name='ActiveEnzymes', index=False)
        pd.read_excel(
            os.path.join('Data', 'proteinAllocationModel_iML1515_EnzymaticData_240730.xlsx'),
            sheet_name='UnusedEnzyme').to_excel(writer, sheet_name='UnusedEnzyme', index=False)
        pd.read_excel(
            os.path.join('Data', 'proteinAllocationModel_iML1515_EnzymaticData_240730.xlsx'),
            sheet_name='Translational').to_excel(writer, sheet_name='Translational', index=False)
        pd.read_excel(
            os.path.join('Data', 'mcPAM_iML1515_EnzymaticData.xlsx'),
            sheet_name='Membrane').to_excel(writer, sheet_name='Membrane', index=False)



    ecoli_pam = setup_ecoli_pam_uniprot(pam_info_file= os.path.join('Data', 'mcPAM_iML1515_EnzymaticData_new.xlsx'))

    ecoli_pam.GLUCOSE_EXCHANGE_RXNID = 'EX_glc__D_e'

    validation_data = set_up_validation_data(c_sources)
    hyperparameters = set_up_hyperparameter(processes, gene_flow_events, filename_extension, num_kcats_to_mutate)

    return PAMParametrizer(pamodel=ecoli_pam,
                     validation_data=validation_data,
                     hyperparameters=hyperparameters,
                     substrate_uptake_id=config.GLUCOSE_EXCHANGE_RXNID,
                     max_substrate_uptake_rate=max_substrate_uptake_rate,
                     min_substrate_uptake_rate=min_substrate_uptake_rate)

if __name__ == "__main__":
    pam_parametrizer = set_up_pamparametrizer(MIN_SUBSTRATE_UPTAKE_RATE, MAX_SUBSTRATE_UPTAKE_RATE,
                         c_sources = ['Glucose'])

    pam_parametrizer.run(remove_subruns=True, binned = 'False')
# for running:
# python -m Scripts.Testing.pam_parametrizer_iML1515
