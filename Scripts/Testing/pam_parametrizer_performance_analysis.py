import argparse
import os
import pandas as pd
from typing import Union
from Scripts.build_pam_parametrizer import run_pam_parametrizer_toymodel
from datetime import date

PARAM_RESULTS_FILE_PATH = os.path.join('Results', 'pam_parametrizer_diagnostics.xlsx')
RESULT_FILE_PATH = os.path.join('Results', f'pam_parametrizer_statistics_{date.today()}.xlsx')

def parse_arguments():
    parser = argparse.ArgumentParser("pam_parametrizer_performance")
    parser.add_argument("--iterations", help="Number of parametrization runs to perform for calculating mean error", type=int)
    args = parser.parse_args()
    return args


def save_pam_parametrizer_results_to_df(iter: int, bin_config: str,
                                        best_individual_df: pd.DataFrame,
                                        computational_time_df: pd.DataFrame):
    best_indiv = pd.read_excel(PARAM_RESULTS_FILE_PATH, sheet_name='Best_Individuals')
    comp_time = pd.read_excel(PARAM_RESULTS_FILE_PATH, sheet_name='Computational_Time')
    final_errors = pd.read_excel(PARAM_RESULTS_FILE_PATH, sheet_name='Final_Errors')

    best_indiv = best_indiv.assign(binned=bin_config, iteration=iter)
    best_indiv = pd.merge(best_indiv, final_errors, on='run_id', how='left')

    print(best_indiv)
    best_individual_df = pd.concat([best_individual_df, best_indiv], axis=0)
    print(best_individual_df)

    comp_time = comp_time.assign(binned=bin_config, iteration=iter)
    computational_time_df = pd.concat([computational_time_df, comp_time], axis=0)
    print(computational_time_df)

    return best_individual_df, computational_time_df

def get_statistics_from_df(df: pd.DataFrame, group_by:Union[list, str],
                           columns:list):
    if isinstance(group_by, str): group_by = [group_by]
    grouped = df.groupby(group_by)

    # Calculate mean, median, and standard deviation for 'ga_error' and 'final_error'
    result = grouped.agg({column: ['mean', 'median', 'std'] for column in columns}).reset_index()
    new_column_names = group_by
    for column in columns:
        new_column_names += [column+'_'+stat for stat in ['mean', 'median', 'std']]

    # Rename the columns for clarity
    result.columns = new_column_names

    return result

if __name__ == '__main__':
    args = parse_arguments()

    if args.iterations is None:
        iterations = 10
    else:
        iterations = args.iterations

    # 0. Initialize result dataframes
    best_individual_df = pd.DataFrame(columns = ['iteration', 'binned',
                                                 'run_id','enzyme_id', 'rxn_id', 'kcat[s-1]',
                                                 'ga_error', 'r_squared'])
    computational_performance_df = pd.DataFrame(columns = ['iteration', 'binned', 'run_id', 'time_s', 'time_h'])

    # 1. Run different configuations of the toy model parametrization for n iterations and save the results
    for iteration in range(iterations):
        print('\n\n###################################################################################')
        print('starting with iteration number ', iteration, ' out of ', iterations, ' iterations\n')
        for configuation in ['False', 'all', 'before']:
            print('Configuration of the parametrizations workflow: ', configuation)
            print('------------------------------------------------------------------------------------------------')
            parametrizer = run_pam_parametrizer_toymodel(binned = configuation)
            best_individual_df, computational_performance_df = save_pam_parametrizer_results_to_df(iteration, configuation,
                                                                               best_individual_df,
                                                                               computational_performance_df)
    # 2. Calculate mean error and stdev for each method and for a single iteration
    error_per_iteration_config = get_statistics_from_df(best_individual_df,
                                                            group_by = ['iteration', 'binned'],
                                                            columns = ['ga_error', 'r_squared'])
    error_per_config = get_statistics_from_df(best_individual_df,
                                                  group_by = ['binned'],
                                                  columns = ['ga_error', 'r_squared'])
    computational_performance_per_config = get_statistics_from_df(computational_performance_df,
                                                                      group_by=['binned'],
                                                                      columns=['time_s'])
    # 3. save the results to excel
    with pd.ExcelWriter(RESULT_FILE_PATH) as writer:
        # Write each DataFrame to a specific sheet
        best_individual_df.to_excel(writer, sheet_name='Best_Individuals', index=False)
        computational_performance_df.to_excel(writer, sheet_name='Computational_Time',
                                                                 index=False)
        error_per_config.to_excel(writer, sheet_name='stats_error_per_config', index=False)
        error_per_iteration_config.to_excel(writer, sheet_name='stats_error_per_iteration', index=False)
        computational_performance_per_config.to_excel(writer, sheet_name='stats_computation_time', index=False)








