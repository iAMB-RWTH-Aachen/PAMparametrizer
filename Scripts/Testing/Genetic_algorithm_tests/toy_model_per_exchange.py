from toy_model import init_toy_parametrization_ga
from PAM_Parametrization.Modules.utils.genetic_algorithm_analysis import get_kcat_error_from_ga

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


def plot_kcats_errors_per_exchange(result_dict:dict,save:bool=False):
    fig = plt.figure(layout="constrained")

    gs0 = GridSpec(1,2,figure=fig)
    gs_error = gs0[0,0]
    gs_kcats = gs0[0,1]

    rxns = list(result_dict.keys())
    errors = [rxn_result['error'] for rxn_result in result_dict.values()]
    kcats = [rxn_result['kcats'] for rxn_result in result_dict.values()]

    ax_error = fig.add_subplot(gs_error)
    ax_error.bar(rxns, errors)
    ax_error.set_title('Error per exchange')
    ax_error.set_xlabel('Reaction')
    ax_error.set_ylabel('error')


    kcat1, kcat2, kcat3 = map(list, zip(*kcats))
    bar_width = 0.15
    index = np.arange(len(rxns))
    ax_kcat = fig.add_subplot(gs_kcats)
    ax_kcat.bar(index, kcat1, color = 'darkred', label = 'E3', width = bar_width)
    ax_kcat.bar(index+bar_width, kcat2, color='darkblue', label='E4', width = bar_width)
    ax_kcat.bar(index+bar_width*2, kcat3, color='purple', label='E5', width = bar_width)
    plt.legend()

    ax_kcat.set_title('Evolution of kcat values')
    ax_kcat.set_xlabel('Reaction')
    ax_kcat.set_ylabel('kcat value [$h^{-1}$]')

    ax_kcat.set_xticks(index + bar_width / 2, rxns)
    fig.set_figwidth(10)
    fig.set_figheight(7)

    if save: plt.savefig('Results/per_exchange.png')
    else: plt.show()

def run_ga_per_exchange(valid_data_per_exchange:pd.DataFrame):
    results = {}
    for rxn, df in valid_data_per_exchange.items():
        print('-----------------------------------------------------------------------------------')
        print('Optimizing to match the following reaction: ', rxn)
        filename = "toy_model_parametrization" + rxn
        ga = init_toy_parametrization_ga(valid_data_df=df,
                                         substrate_uptake_rates=df['R1_ub'].to_list(),
                                         filename_save=filename
                                         )
        ga.start()
        result_dict = get_kcat_error_from_ga(ga)

        print('kcats: ', result_dict['kcats'], '\nerror: ', result_dict['error'])
        print('-----------------------------------------------------------------------------------')

        results  = {**results, rxn: result_dict}

    return results

def restart_ga(rxns:list):
    print('-----------------------------------------------------------------------------------')
    print('Restarting genetic algorithm with populations optimized for the following reactions', ' '.join(rxns))
    filenames = ["Results/toy_model_parametrization" + rxn + '.json' for rxn in rxns]
    final_results_file_path = "toy_model_parametrization_restart_per_exchange"
    ga = init_toy_parametrization_ga(filename_save=final_results_file_path)
    ga.restart(filenames)
    result_dict = get_kcat_error_from_ga(ga)
    print('kcats: ', result_dict['kcats'], '\nerror: ', result_dict['error'])
    print('-----------------------------------------------------------------------------------')

if __name__ == '__main__':
    # start kcat: [1, 0.5, 1, 0.5 ,0.45, 1.5]
    # aimed end kcat: [1, 0.5, 5, 0.1, 0.25, 1.5]

    # get the validation data file
    DATA_DIR = os.path.join(os.path.split(os.getcwd())[0], 'Data')
    RESULT_DF_FILE = os.path.join(DATA_DIR, 'toy_model_simulations_ga.csv')
    valid_data_df = pd.read_csv(RESULT_DF_FILE)

    # separate per exchange
    valid_data_per_exchange = {rxn: valid_data_df[['R1_ub', rxn]] for rxn in valid_data_df.columns[2:]}
    results = run_ga_per_exchange(valid_data_per_exchange)
    #to do restart algorithm using combination of the populations of the previous run
    #to do try uniform sampling
    plot_kcats_errors_per_exchange(result_dict=results, save=True)

    restart_ga(list(valid_data_per_exchange.keys()))