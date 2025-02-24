import os.path
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import pandas as pd
from Scripts.i3_analysis.PAMparametrizer_progress_cleaned_figure import recreate_progress_plot
from Scripts.i2_parametrization.mcpam import set_up_pamparametrizer as set_up_pamparametrizer_mcpam
from Scripts.i2_parametrization.pam_parametrizer_iML1515 import set_up_pamparametrizer as set_up_pamparametrizer_ecoli



result_file_names = ['pam_parametrizer_statistics_2024-08-14.xlsx']

def create_progress_plots(result_file_names,
                          output_file_path,
                          recreate_plot=False,
                          pam_parametrizer_setup = set_up_pamparametrizer_ecoli):
    max_iteration = 10
    fig, axs = plt.subplots(ncols=len(result_file_names), sharey='row')
    with pd.ExcelWriter(output_file_path, engine='openpyxl') as writer:
        pd.DataFrame(columns=['something']).to_excel(writer)

    for i,file in enumerate(result_file_names):
        best_indiv_all_rounds = pd.read_excel(os.path.join('Results', '2_parametrization', 'diagnostics', file),
                                              sheet_name='Best_Individuals')
        print(best_indiv_all_rounds)
        best_indiv_all_rounds = best_indiv_all_rounds.drop_duplicates(subset=['enzyme_id', 'rxn_id', 'ga_error', 'direction'],
                                                                      keep = 'first')
        best_indiv_grouped = best_indiv_all_rounds.groupby('iteration')
        config = best_indiv_all_rounds.iloc[0]['binned']

        for iteration, best_indiv_df in best_indiv_grouped:
            cmap = plt.get_cmap('viridis')
            color = to_hex(cmap(iteration / (max_iteration + 1)))

            if recreate_plot:
                fig_file_path = os.path.join('Scripts', 'Results', 'i3_analysis',
                                             f'pam_parametrizer_progess_cleaned_{config}_{iteration}_mcpam.png')
                error_df = recreate_progress_plot(best_indiv_df, fig_file_path,
                                                  return_error_df = True,
                                                  setup_pam_parametrizer = pam_parametrizer_setup)
                with pd.ExcelWriter(output_file_path, engine = 'openpyxl', mode='a', if_sheet_exists='replace') as writer:
                    error_df.to_excel(writer, sheet_name=f'error_{config}_{iteration}')

                if len(result_file_names)==1:
                    axs.plot(error_df.run_id, error_df.error, color=color)
                else:
                    axs[i].plot(error_df.run_id, error_df.error, color = color)


            else:
                axs[i].plot(best_indiv_df.run_id, best_indiv_df.ga_error, color=color)

    plt.show()

def main_mcecoli():
    result_file_names = ['pam_parametrizer_diagnostics_mciML1515_1.xlsx']

    output_file_path = os.path.join('Results','3_analysis', 'error_dataframe_for_checking_MCPAM.xlsx')

    #
    fig_file_path = os.path.join('Results', '3_analysis', 'pam_parametrizer_progess_cleaned_mcpam.png')

    create_progress_plots(result_file_names,
                           output_file_path,
                           recreate_plot = True,
                           pam_parametrizer_setup=set_up_pamparametrizer_mcpam)

def main_ecoli():
    output_file_path = os.path.join('Results','3_analysis', 'error_dataframe_for_checking.xlsx')
    create_progress_plots(result_file_names,output_file_path, True)

if __name__ == '__main__':
    main_mcecoli()
