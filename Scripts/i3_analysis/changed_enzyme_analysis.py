import os.path
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import pandas as pd

result_file_names = ['pam_parametrizer_statistics_2024-06-14.xlsx']

def visualize_changed_enzymes(result_file_names, output_file_path, recreate_plot=False):
    max_iteration = 0
    fig, axs = plt.subplots(ncols=len(result_file_names), sharey='row')
    best_indiv_all_rounds = pd.read_excel(os.path.join('Results', file), sheet_name='Best_Individuals')
    best_indiv_all_rounds = best_indiv_all_rounds.drop_duplicates(
        subset=['enzyme_id', 'rxn_id', 'r_squared', 'direction'], keep='first')
    best_indiv_grouped = best_indiv_all_rounds.groupby('iteration')
    config = best_indiv_all_rounds.iloc[0]['binned']


    for iteration, best_indiv_df in best_indiv_grouped:
        cmap = plt.get_cmap('viridis')
        color = to_hex(cmap(iteration / (max_iteration + 1)))
