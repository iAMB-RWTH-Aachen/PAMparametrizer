import os
from datetime import date
import subprocess
import pandas as pd

from Scripts.i2_parametrization.pam_parametrizer_toy_model import set_up_pamparametrizer
from Scripts.i2_parametrization.pam_parametrizer_performance_analysis import run_parametrization_workflow

RESULT_FILE_PATH = os.path.join('Results', 'i2_parametrization', f'pam_parametrizer_statistics_{date.today()}.xlsx')



def test_pam_performance_analysis_runs_from_terminal_command():
    # Arrange
    command = 'python3 -m Scripts.i2_parametrization.pam_parametrizer_performance_analysis --hyper_processes 1 --iterations 2 --model toy'

    RESULT_PARAM_FILE_PATH_XLS = os.path.join('Results', 'i2_parametrization', 'diagnostics',
                                              f'pam_parametrizer_diagnostics_0.xlsx')
    RESULT_PARAM_FILE_PATH_PNG = os.path.join('Results', 'i2_parametrization', 'progress',
                                              f'pam_parametrizer_progress_0.png')

    # Apply
    subprocess.run(command, shell=True)
    # Assert
    assert os.path.exists(RESULT_FILE_PATH)
    [os.remove(path) for path in [RESULT_FILE_PATH, RESULT_PARAM_FILE_PATH_PNG]]#, RESULT_PARAM_FILE_PATH_XLS]]

def test_pam_performance_analysis_saves_results_correctly():
    # Arrange
    best_individual_df = pd.DataFrame(columns=['iteration', 'binned',
                                               'run_id', 'enzyme_id', 'rxn_id', 'kcat[s-1]',
                                               'ga_error', 'r_squared'])
    computational_performance_df = pd.DataFrame(columns=['iteration', 'binned', 'run_id', 'time_s', 'time_h'])
    configurations = ['False', 'before', 'all']
    iterations = 2
    processes = 1
    gene_flow_events = 1
    num_kcats_to_mutate= 3
    # files_to_remove = [os.path.join('Results',f'pam_parametrizer_diagnostics_{config}.xlsx') for config in configurations]
    files_to_remove = [os.path.join('Results','i2_parametrization', 'progress',
                                    f'pam_parametrizer_progress_{config}.png') for config in configurations]

    # Apply
    for iteration in range(iterations):
        best_individual_df, computational_performance_df = run_parametrization_workflow(iteration, iterations,
                                 configurations, set_up_pamparametrizer,processes, gene_flow_events,
                                 num_kcats_to_mutate,
                                 best_individual_df, computational_performance_df, 1e-3, 0.1)

    # best_individual_df_group1 = best_individual_df.groupby(['iteration', 'run_id', 'binned'])
    best_individual_df_group1 =best_individual_df[(best_individual_df.run_id == 1) & (best_individual_df.iteration == 1)
                                                  & (best_individual_df.binned == 'all')]
    best_individual_df_no_f = best_individual_df_group1[best_individual_df_group1.direction != 'f']
    best_individual_df_no_b = best_individual_df_group1[best_individual_df_group1.direction != 'b']

    # Assert
    assert all([len(df) >= num_kcats_to_mutate for df in [best_individual_df_no_b, best_individual_df_no_f]])
    [os.remove(path) for path in files_to_remove]
