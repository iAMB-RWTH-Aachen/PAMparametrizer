import os
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

from Modules.PAM_parametrizer import ValidationData, HyperParameters, ParametrizationResults, SectorConfig
from Modules.PAM_parametrizer import PAMParametrizer
from Scripts.pam_generation import setup_toy_pam

MAX_SUBSTRATE_UPTAKE_RATE = 0.1
MIN_SUBSTRATE_UPTAKE_RATE = 0.001

def set_up_validation_data():
    DATA_DIR = os.path.join('tests', 'data')
    RESULT_DF_FILE = os.path.join(DATA_DIR, 'toy_model_simulations_ga.csv')
    valid_data_df = pd.read_csv(RESULT_DF_FILE)

    validation_data = ValidationData(valid_data_df, 'R1', [MIN_SUBSTRATE_UPTAKE_RATE,MAX_SUBSTRATE_UPTAKE_RATE])
    validation_data.sector_configs = {'TranslationalProteinSector':SectorConfig(
            sectorname = 'TranslationalProteinSector',
            slope = 0.01*1e-3,
            intercept = 0.01*1e-3,
            substrate_range = [-1e-3,-2*1e-3]
        )}
    validation_data._reactions_to_plot = ['R1', 'R7', 'R8', 'R9']
    validation_data._reactions_to_validate = ['R1', 'R7', 'R8', 'R9']
    return validation_data

def set_up_hyperparameter(processes: int,
                          gene_flow_events:int,
                          filename_extension:str = 'toy_',
                          num_kcats_to_mutate:int=3):
    hyperparams = HyperParameters
    hyperparams.threshold_iteration = 4
    hyperparams.threshold_error = 0.95
    hyperparams.number_of_kcats_to_mutate = num_kcats_to_mutate
    hyperparams.filename_extension = filename_extension

    hyperparams.genetic_algorithm_hyperparams['number_generations'] = 2
    hyperparams.genetic_algorithm_hyperparams['number_gene_flow_events'] = gene_flow_events
    hyperparams.genetic_algorithm_hyperparams['processes'] = processes

    hyperparams.genetic_algorithm_filename_base = 'genetic_algorithm_run_toy_'
    hyperparams.genetic_algorithm_hyperparams['print_progress'] = False
    return hyperparams

def set_up_toy_model(kcat_fwd:list = [1, 0.5, 1, 0.5, 0.45, 1.5]):
    # kcat_fwd = [1, 0.5, 5, 0.1, 0.25, 1.5] #the 'final' dataset
    toy_pam = setup_toy_pam(kcat_fwd=kcat_fwd)
    toy_pam.name = 'toy_model'
    return toy_pam

def run_simulations(pamodel, substrate_rates):
    result_df = pd.DataFrame(columns= ['R1_ub','R1', 'R7', 'R8', 'R9'])

    for substrate in substrate_rates:
        pamodel.change_reaction_bounds(rxn_id='R1',
                                       lower_bound=0, upper_bound=substrate)
        print('Running simulations with ', substrate, 'mmol/g_cdw/h of substrate going into the system')
        pamodel.optimize()
        if pamodel.solver.status == 'optimal' and pamodel.objective.value>0:
            results_row = []
            for rxn in ['R1', 'R7', 'R8', 'R9']:
                results_row += [pamodel.reactions.get_by_id(rxn).flux]

            result_df.loc[len(result_df)] = [substrate] + results_row
    return result_df

def set_up_pamparametrizer(min_substrate_uptake_rate:float, max_substrate_uptake_rate: float,
                           pam_info_file = None,
                           kcat_increase_factor = 1,
                           processes: int =2,
                           gene_flow_events: int = 2,
                           filename_extension:str = 'toy',
                           num_kcats_to_mutate:int =3,
                           kcat_fwd: list = [1, 0.5, 1, 0.5, 0.45, 1.5], sensitivity:bool = True, c_sources=None):

    params = {'pamodel': set_up_toy_model(kcat_fwd),
              'validation_data': set_up_validation_data(),
              'hyperparameters': set_up_hyperparameter(processes, gene_flow_events, filename_extension, num_kcats_to_mutate),
              'substrate_uptake_id': 'R1',
              'max_substrate_uptake_rate': max_substrate_uptake_rate,
              'min_substrate_uptake_rate': min_substrate_uptake_rate,
                'sensitivity': sensitivity}
    if not sensitivity:
        params['enzymes_to_evaluate'] = {'E3':{'reaction':'R3','kcat':1, 'sensitivity':0.5}, #should become 5
                           'E4':{'reaction':'R4','kcat':0.5, 'sensitivity':0.2},#should become 0.1
                           'E5':{'reaction':'R5','kcat':0.45, 'sensitivity':0.1}},#should become 0.25

    return PAMParametrizer(**params)

if __name__ == "__main__":
    # pam_parametrizer = set_up_pamparametrizer('False')
    pam_parametrizer = set_up_pamparametrizer(MIN_SUBSTRATE_UPTAKE_RATE, MAX_SUBSTRATE_UPTAKE_RATE)
    pam_parametrizer.run(binned = 'before')
    # pam_parametrizer_no_senz = set_up_pamparametrizer(MIN_SUBSTRATE_UPTAKE_RATE, MAX_SUBSTRATE_UPTAKE_RATE,
    #                                                   sensitivity= False)
    # pam_parametrizer_no_senz.run()
# for running:
# python -m Scripts.i2_parametrization.pam_parametrizer_toy_model
