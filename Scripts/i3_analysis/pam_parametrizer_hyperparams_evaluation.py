import pandas as pd
import time
import os

from Scripts.Testing.pam_parametrizer_iML1515 import set_up_pamparametrizer as set_up_pamparametrizer_iml

def scan_number_enzymes_to_mutate(number_of_iterations:int = 3,
                                  number_enzymes_list:list = [5, 10, 15, 20, 100],
                                  result_file_path: str = os.path.join('Results', 'num_enzymes_to_mutate_scan_iml1515.xlsx')) -> pd.DataFrame:
    num_enzymes_error_time_df = pd.DataFrame(columns=['number_enzymes', 'final_error', 'comp_time', 'iteration'])
    for num_enzymes in number_enzymes_list:
        print('\n-------------------------------------------------------------------------------------------')
        print(f'Running PAMparametrizer while mutating {num_enzymes} kcats in a single iteration')
        for i in range(1,number_of_iterations+1):
            print(f'==========================\nRunning replicate {i} out of {number_of_iterations}\n\n')
            pam_parametrizer = set_up_pamparametrizer_iml(min_substrate_uptake_rate=-11,
                                                          max_substrate_uptake_rate=-1,
                                                          processes=2, gene_flow_events=2,
                                                          filename_extension=f'iML1515_num_enzyme_{num_enzymes}_{i}',
                                                          num_kcats_to_mutate=num_enzymes,
                                                          threshold_iteration=3,
                                                          kcat_increase_factor=3)

            # need to reset best individual and computational performance df
            pam_parametrizer.parametrization_results.best_individuals = pd.DataFrame(
                columns=['run_id', 'enzyme_id', 'direction', 'rxn_id', 'kcat[s-1]', 'ga_error'])
            pam_parametrizer.parametrization_results.computational_time = pd.DataFrame(
                columns=['run_id', 'time_s', 'time_h'])

            # reset final errors for correct saving
            pam_parametrizer.parametrization_results.final_errors = pd.DataFrame(columns=['run_id', 'r_squared'])

            start = time.perf_counter()
            pam_parametrizer.run(remove_subruns=True, binned='False')
            total_runtime = time.perf_counter() - start

            num_enzymes_error_time_df.loc[len(num_enzymes_error_time_df)] = [num_enzymes, pam_parametrizer.final_error, total_runtime, i]
        save_intermediate_result_df(num_enzymes_error_time_df, result_file_path)

def save_intermediate_result_df(df:pd.DataFrame, file_path:str) -> None:
    write_mode = 'w'
    kwargs = {}

    if os.path.isfile(file_path):
        write_mode = 'a'
        kwargs = {'if_sheet_exists': 'replace'}

    with pd.ExcelWriter(file_path, mode=write_mode, engine='openpyxl', **kwargs) as writer:
        df.to_excel(writer, index=False)

if __name__ == '__main__':
    scan_number_enzymes_to_mutate(number_enzymes_list=[100])

