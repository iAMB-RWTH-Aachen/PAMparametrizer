import argparse
import os
import pandas as pd
from typing import Union
from Scripts.i2_parametrization.pam_parametrizer_iML1515 import set_up_pamparametrizer as set_up_pamparametrizer_iml
from Scripts.i2_parametrization.pam_parametrizer_toy_model import set_up_pamparametrizer as set_up_pamparametrizer_toy
from Scripts.i2_parametrization.pam_parametrizer_ecolicore import set_up_pamparametrizer as set_up_pamparametrizer_core
from Modules.utils.pam_generation import set_up_pam


from datetime import date

RESULT_FILE_PATH = os.path.join('Results', '3_parametrization',  f'pam_parametrizer_statistics_{date.today()}.xlsx')

def parse_arguments():
    parser = argparse.ArgumentParser("pam_parametrizer_performance")
    parser.add_argument("--model",
                        help="which model to use can be 'ecolicore', 'ecoli' or 'toy'",
                        type=str)
    parser.add_argument("--pam_info_file",
                        help="path to the file containing information about the parameters to build the pam",
                        type=str)
    parser.add_argument("--configuration",
                        help="Binning configuration of the parameterization workflow can be 'all', 'before' or 'False'",
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


def save_pam_parametrizer_results_to_df(iter: int, bin_config: str,
                                        best_individual_df: pd.DataFrame,
                                        computational_time_df: pd.DataFrame,
                                        result_file_path:str):
    best_indiv = pd.read_excel(result_file_path, sheet_name='Best_Individuals')
    comp_time = pd.read_excel(result_file_path, sheet_name='Computational_Time')
    final_errors = pd.read_excel(result_file_path, sheet_name='Final_Errors')

    best_indiv = best_indiv.assign(binned=bin_config, iteration=iter)
    best_indiv = pd.merge(best_indiv, final_errors, on='run_id', how='left')

    best_individual_df = pd.concat([best_individual_df, best_indiv], axis=0).reset_index(drop=True)#.drop_duplicates(
        #subset=['enzyme_id', 'rxn_id', 'r_squared', 'direction'], keep = 'first')

    comp_time = comp_time.assign(binned=bin_config, iteration=iter)
    computational_time_df = pd.concat([computational_time_df, comp_time], axis=0)

    # os.remove(result_file_path)
    return best_individual_df, computational_time_df

def get_statistics_from_df(df: pd.DataFrame, group_by:Union[list, str],
                           columns:list):
    if isinstance(group_by, str): group_by = [group_by]
    grouped = df.groupby(group_by)

    # Calculate mean, median, and standard deviation for 'ga_error' and 'final_error'
    result = grouped.agg({column: ['mean', 'median', 'std', 'min', 'max'] for column in columns}).reset_index()
    new_column_names = group_by
    for column in columns:
        new_column_names += [column+'_'+stat for stat in ['mean', 'median', 'std', 'min', 'max']]

    # Rename the columns for clarity
    result.columns = new_column_names

    return result

def run_parametrization_workflow(iteration, iterations,
                                 configurations, set_up_pamparametrizer,
                                 processes,
                                 gene_flow_events, num_kcats_to_mutate,
                                 best_individual_df, computational_performance_df,
                                 pam_info_file: str,
                                 min_substrate_uptake = -11, max_substrate_uptake = -0.1):
    print('\n\n###################################################################################')
    print('starting with iteration number ', iteration, ' out of ', iterations, ' iterations\n')
    for configuration in configurations:
        print('Working on iteration number', iteration, 'out of ',iterations)
        print('Configuration of the parametrizations workflow: ', configuration)
        print('------------------------------------------------------------------------------------------------')
        parametrizer = set_up_pamparametrizer(min_substrate_uptake, max_substrate_uptake,
                                              pam_info_file = pam_info_file,
                                              processes=processes,
                                              gene_flow_events=gene_flow_events,
                                              filename_extension= str(iteration)+'_1',
                                              num_kcats_to_mutate = num_kcats_to_mutate,
                                              c_sources = ['Glucose'],
                                              kcat_increase_factor=3#, 'Glycerol', 'Acetate']
                                              #['Glycerol', 'Glucose', 'Acetate', 'Pyruvate', 'Gluconate', 'Succinate', 'Galactose', 'Fructose']
                                              )
        pamodel_copy = parametrizer.pamodel.copy(copy_with_pickle = True)
        parametrizer.run(binned=configuration)

        #need to reset best individual and computational performance df
        parametrizer.parametrization_results.best_individuals = pd.DataFrame(columns=['run_id', 'enzyme_id', 'direction', 'rxn_id', 'kcat[s-1]', 'ga_error'])
        parametrizer.parametrization_results.computational_time = pd.DataFrame(columns=['run_id', 'time_s', 'time_h'])

        #reset final errors for correct saving
        parametrizer.parametrization_results.final_errors = pd.DataFrame(columns=['run_id', 'r_squared'])
        if len(configurations)>1:
            parametrizer.pamodel = pamodel_copy

        results_file_path = os.path.join('Results', 'i2_parametrization', 'diagnostics',
                                         f'pam_parametrizer_diagnostics_{str(iteration)}.xlsx')
        best_individual_df, computational_performance_df =  save_pam_parametrizer_results_to_df(iteration,configuration,
                                                   best_individual_df,computational_performance_df,
                                                   results_file_path)
    return best_individual_df, computational_performance_df

def analyse_parametrizer_performance():
    min_substrate = -11
    max_substrate = -0.1

    args = parse_arguments()
    pam_info_file =  args.pam_info_file
    iterations = args.iterations
    configuration = args.configuration
    processes = args.hyper_processes
    gene_flow_events = args.hyper_gfe

    if iterations is None:
        iterations = 5
    if pam_info_file is None:
        pam_info_file = os.path.join('Results','1_preprocessing', 'proteinAllocationModel_iML1515_EnzymaticData_250214.xlsx')
    if configuration is None:
        configurations = ['False', 'all', 'before']
    else:
        configurations = [configuration]
    if processes is None:
        processes = 2
    if gene_flow_events is None:
        gene_flow_events = processes
    if args.model == 'toy':
        set_up_pamparametrizer = set_up_pamparametrizer_toy
        min_substrate = 1e-3
        max_substrate = 0.1
    elif args.model == 'ecolicore':
        set_up_pamparametrizer = set_up_pamparametrizer_core
    else:
        set_up_pamparametrizer = set_up_pamparametrizer_iml

    # 0. Initialize result dataframes
    best_individual_df = pd.DataFrame(columns=['iteration', 'binned',
                                               'run_id', 'enzyme_id', 'rxn_id', 'kcat[s-1]',
                                               'ga_error', 'r_squared'])
    computational_performance_df = pd.DataFrame(columns=['iteration', 'binned', 'run_id', 'time_s', 'time_h'])

    # 1. Run different configuations of the toy model parametrization for n iterations and save the results
    for iteration in range(iterations):
        best_individual_df, computational_performance_df = run_parametrization_workflow(iteration+1, iterations,
                                 configurations, set_up_pamparametrizer,
                                 processes,
                                 gene_flow_events, 10,
                                 best_individual_df, computational_performance_df,
                                pam_info_file, min_substrate_uptake=min_substrate,
                                 max_substrate_uptake=max_substrate)

    # 2. Calculate mean error and stdev for each method and for a single iteration
    error_per_iteration_config = get_statistics_from_df(best_individual_df,
                                                        group_by=['iteration', 'binned'],
                                                        columns=['ga_error', 'r_squared'])
    # 3. save the results to excel
    with pd.ExcelWriter(RESULT_FILE_PATH) as writer:
        # Write each DataFrame to a specific sheet
        best_individual_df.to_excel(writer, sheet_name='Best_Individuals', index=False)
        computational_performance_df.to_excel(writer, sheet_name='Computational_Time',
                                              index=False)
        # error_per_config.to_excel(writer, sheet_name='stats_error_per_config', index=False)
        error_per_iteration_config.to_excel(writer, sheet_name='stats_error_per_iteration', index=False)
        # computational_performance_per_config.to_excel(writer, sheet_name='stats_computation_time', index=False)


if __name__ == '__main__':
    analyse_parametrizer_performance()







