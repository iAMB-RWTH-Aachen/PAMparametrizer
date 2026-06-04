import argparse
import os
import pandas as pd
from collections.abc import Iterable

from Scripts.i2_parametrization.pam_parametrizer_iML1515 import set_up_pamparametrizer, set_up_validation_data

from PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_uniform import FitnessEvaluation
from PAMparametrizer.utils.pamparametrizer_analysis import calulate_pearson_correlation_simulation_vs_experiment


class FitnessEvaluationPearson(FitnessEvaluation):
    def _calculate_simulation_error(self, flux_df: pd.DataFrame, substrate_reaction:str):
        correlation, p_value = calulate_pearson_correlation_simulation_vs_experiment(
            validation_df=self.valid_data[substrate_reaction],
            flux_df=flux_df,
            rxns_to_validate=self.reactions_with_data[substrate_reaction],
            substr_rxn=substrate_reaction+'_ub',
            substrate_sim='substrate',

        )

        return correlation

def parse_arguments():
    parser = argparse.ArgumentParser("pam_parametrizer_performance")
    parser.add_argument("--pam_info_file",
                        help="path to the file containing information about the parameters to build the pam",
                        type=str)
    parser.add_argument("--iterations",
                        help="Number of parametrization runs to perform for calculating mean error",
                        type=int)
    parser.add_argument("--hyper_processes",
                        help="Number of parallel workers for parametrization workflow",
                        type=int)
    parser.add_argument("--hyper_gfe",
                        help="Number of gene flow events, i.e. merging of multiple populations independently evolved on parallel workers",
                        type=int)

    args = parser.parse_args()
    return args


def run_parametrization_workflow(iteration, iterations,set_up_pamparametrizer,
                                 processes,
                                 gene_flow_events, num_kcats_to_mutate,
                                 pam_info_file: str,
                                 min_substrate_uptake = -11, max_substrate_uptake = -0.1):
    print('\n\n###################################################################################')
    print('starting with iteration number ', iteration, ' out of ', iterations, ' iterations\n')
    parametrizer = set_up_pamparametrizer(min_substrate_uptake, max_substrate_uptake,
                                              pam_info_file = pam_info_file,
                                              processes=processes,
                                              gene_flow_events=gene_flow_events,
                                              filename_extension= str(iteration),
                                              num_kcats_to_mutate = num_kcats_to_mutate,
                                              c_sources = ['Glucose'],
                                              threshold_iteration = 10,
                                              kcat_increase_factor=7#, 'Glycerol', 'Acetate']
                                              #['Glycerol', 'Glucose', 'Acetate', 'Pyruvate', 'Gluconate', 'Succinate', 'Galactose', 'Fructose']
                                              )
    parametrizer.hyperparameters.genetic_algorithm_hyperparams['fitness_class'] = FitnessEvaluationPearson
    parametrizer.run()



def analyse_parametrizer_performance():
    min_substrate = -11
    max_substrate = -0.1

    args = parse_arguments()
    pam_info_file =  args.pam_info_file
    iterations = args.iterations
    processes = args.hyper_processes
    gene_flow_events = args.hyper_gfe

    if iterations is None:
        iterations = 5
    if pam_info_file is None:
        pam_info_file = os.path.join('Results','1_preprocessing',
                                     'proteinAllocationModel_iML1515_EnzymaticData_250912.xlsx')
    if processes is None:
        processes = 2
    if gene_flow_events is None:
        gene_flow_events = processes

    for iteration in range(iterations):
        run_parametrization_workflow(
            iteration+1, iterations,set_up_pamparametrizer,
            processes,
            gene_flow_events, 40,
            pam_info_file, min_substrate_uptake=min_substrate,
            max_substrate_uptake=max_substrate)

if __name__ == '__main__':
    analyse_parametrizer_performance()
    # vd = set_up_validation_data(['Glucose'], 'Results/1_preprocessing/proteinAllocationModel_iML1515_EnzymaticData_250912.xlsx')
    # print(vd)






