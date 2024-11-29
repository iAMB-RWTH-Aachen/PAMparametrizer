import os
import pandas as pd
import warnings

warnings.filterwarnings("ignore")

from PAModelpy.configuration import Config

from Modules.PAM_parametrizer import ValidationData, HyperParameters, ParametrizationResults
from Modules.PAM_parametrizer import PAMParametrizer
from Scripts.pam_generation_uniprot_id import setup_yeast_pam, increase_kcats_in_parameter_file

MAX_SUBSTRATE_UPTAKE_RATE = -0.1
MIN_SUBSTRATE_UPTAKE_RATE = -20
config = Config()
config.reset()

def set_up_validation_data(model,csources: list=None) -> list[ValidationData]:
    condition2uptake = {'Maltose': 'r_1931', 'Glucose': 'r_1714',
                        'Galactose': 'r_1710', 'Trehalose': 'r_1650'}
    if csources is None: csources = list(condition2uptake.keys())

    model_reactions = [rxn.id for rxn in model.reactions]

    VALID_DATA_PATH = os.path.join('Data', 'Scerevisiae_phenotypes', 'scerevisiae_phenotypes.xlsx')
    filtered_condition2uptake = {}
    for csource, uptake_rxn in condition2uptake.items():
        if (uptake_rxn in model_reactions) and (csource in csources):
            filtered_condition2uptake[csource] = uptake_rxn

    # Load the data from the sheets
    exchanges_chemostat = pd.read_excel(VALID_DATA_PATH, 'exchanges_chemostat').drop(['RQ'], axis = 1)
    exchanges_batch = pd.read_excel(VALID_DATA_PATH, 'exchanges_batch').drop(['Substrate', 'RQ'], axis = 1)

    exchanges = pd.concat([exchanges_batch, exchanges_chemostat], axis =0)


    #make validation data objects
    validation_data_objects = []
    for c_uptake_id in filtered_condition2uptake.values():
        if c_uptake_id in exchanges.columns:
            #get all rows with uptake of this carbon source and drop all empty columns
            valid_data_df = exchanges[~exchanges[c_uptake_id].isnull()].dropna(axis=1, how='all')

        valid_data_df[c_uptake_id + '_ub'] = valid_data_df[c_uptake_id]

        validation_data = ValidationData(valid_data_df, c_uptake_id, [-20, -0.1])
        #validate only exchange rates and growth rate
        validation_data._reactions_to_validate = [rxn for rxn in valid_data_df.columns if (rxn[-3:]!="_ub")]
        validation_data._reactions_to_plot = ['r_2111', 'r_1761', 'r_1672', 'r_1992', 'r_2033']


        if c_uptake_id == 'r_1714': #glucose uptake
            validation_data.translational_sector_config = {
                'slope': model.sectors.get_by_id('TranslationalProteinSector').tps_mu[0],
                'intercept': model.sectors.get_by_id('TranslationalProteinSector').tps_0[0]
            }

        validation_data_objects.append(validation_data)

    return validation_data_objects

def set_up_hyperparameter(processes: int,
                          gene_flow_events:int,
                          filename_extension:str,
                          num_kcats_to_mutate:int = 4,
                          threshold_iteration:int=10):
    hyperparams = HyperParameters
    hyperparams.threshold_iteration = threshold_iteration
    hyperparams.number_of_kcats_to_mutate = num_kcats_to_mutate
    hyperparams.filename_extension = filename_extension
    hyperparams.genetic_algorithm_filename_base += filename_extension

    hyperparams.genetic_algorithm_hyperparams['processes'] = processes
    hyperparams.genetic_algorithm_hyperparams['number_gene_flow_events'] = gene_flow_events
    hyperparams.genetic_algorithm_hyperparams['number_generations'] = 5
    hyperparams.genetic_algorithm_filename_base = 'genetic_algorithm_run_yeast_'
    hyperparams.genetic_algorithm_hyperparams['print_progress'] = True

    hyperparams.genetic_algorithm_hyperparams['error_weights'] = {'r_1761':3,
                                                                  config.BIOMASS_REACTION: 7}
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
                           filename_extension:str = 'yeast9',
                           num_kcats_to_mutate: int =10,
                           threshold_iteration: int = 10,
                           c_sources:list = ['Glucose'],
                           kcat_increase_factor: int = 1):

    pam_info_file_path_new = os.path.join(
        'Results', '1_preprocessing', 'proteinAllocationModel_yeast9_EnzymaticData_240903_multi.xlsx')
    increase_kcats_in_parameter_file(kcat_increase_factor,
                                     pam_info_file_path_ori=os.path.join(
                                         'Results', '1_preprocessing', 'proteinAllocationModel_yeast9_EnzymaticData_240903.xlsx'),
                                     pam_info_file_path_out=pam_info_file_path_new)

    yeast_pam = setup_yeast_pam()

    validation_data = set_up_validation_data(yeast_pam, c_sources)
    hyperparameters = set_up_hyperparameter(processes, gene_flow_events,
                                            filename_extension,
                                            num_kcats_to_mutate, threshold_iteration)

    return PAMParametrizer(pamodel=yeast_pam,
                     validation_data=validation_data,
                     hyperparameters=hyperparameters,
                     substrate_uptake_id='r_1714',
                     max_substrate_uptake_rate=max_substrate_uptake_rate,
                     min_substrate_uptake_rate=min_substrate_uptake_rate)

if __name__ == "__main__":
    pam_parametrizer = set_up_pamparametrizer(MIN_SUBSTRATE_UPTAKE_RATE, MAX_SUBSTRATE_UPTAKE_RATE,
                                              threshold_iteration= 5, c_sources = ['Glucose'],
                                              kcat_increase_factor=1, processes=2, gene_flow_events=2)#, 'Succinate', 'Fructose','Octanoate','m-Xylene','Toluene','Benzoate'])
    #
    pam_parametrizer.run(remove_subruns=True, binned = 'False')
# for running:
# python -m Scripts.i2_parametrization.pam_parametrizer_iML1515