import numpy as np
import time
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

#ignoring infeasible warnings for cleaner output
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

os.chdir('../Testing/Genetic_algorithm_tests/')


from PAM_Parametrization.Scripts.Testing.Genetic_algorithm_tests.toy_model import init_toy_parametrization_ga

sys.path.append('../../Visualization')

hyperparameters = {'mutation_probability':np.arange(0,1, 0.1), # probability with which an individual (solution) is mutated in a generation
    'mutation_rate':np.arange(0.0001,0.01,0.005 ), # probability with which an attribute (e.g. gene) of an individual is mutated
    'population_size': np.arange(4,50,2), # number of individuals (solution) per population
    'crossover_probability':np.arange(0,1, 0.1), # probability with which two indivduals/offsprings are crossed over
    'number_generations':np.arange(4,50,2), # number of consecutive generations per gene flow event
    'number_gene_flow_events':np.arange(1,10,1), # number of gene flow events, i.e. merging of multiple
                                 # populations independently evolved on parallel workers
    'init_attribute_probability':np.arange(0.0001,0.01,0.005),# probability with which attributes of initial individuals are mutated],
    'sigma_denominator': np.arange(1,10,1) # determines spread of distribution to sample (defined as kcat/sd_denominator). The larger the value, the smaller the spread
                   }

hyperparameters_defaults = {'mutation_probability' : 0.5, # probability with which an individual (solution) is mutated in a generation
    'mutation_rate' : 0.01, # probability with which an attribute (e.g. gene) of an individual is mutated
    'population_size' : 10, # number of individuals (solution) per population
    'crossover_probability' : 0.8, # probability with which two indivduals/offsprings are crossed over
    'number_generations' : 20, # number of consecutive generations per gene flow event
    'number_gene_flow_events' : 2, # number of gene flow events, i.e. merging of multiple
                                 # populations independently evolved on parallel workers
    'init_attribute_probability':0.001,
    'sigma_denominator':10}

def run_GA_custom_hyperparams(hyperparams:dict = hyperparameters_defaults):
    ga = init_toy_parametrization_ga(**hyperparams)
    ga.start()
    result_path = str(ga.folderpath_save.joinpath(ga.filename_save)) + '.xlsx'
    results = pd.read_excel(result_path, sheet_name='final_population').sort_values(by=['fitness_weighted_sum'],ascending=False)

    best_indiv_kcats = [float(kcat) for kcat in results.iloc[0][['attributes']].iloc[0].split(',')]
    best_r2 = results.iloc[0][['fitness_weighted_sum']].iloc[0]

    return {'kcats': best_indiv_kcats, 'error':best_r2}

def hyperparam_range(hyperparameter: str, range: np.arange, save:bool=False):
    results = {}
    for param in range:
        #converting np.int or np.float to normal int or float
        if 'number' in hyperparameter or 'size' in hyperparameter:
            param = int(param)
        else:
            param = float(param)
        #print progress
        print('-------------------------------------------------------------------------------------')
        print('Changing ', hyperparameter, 'to ', param)
        print('-------------------------------------------------------------------------------------')

        #set the parameter dict with the parameter to test
        param_dict = hyperparameters_defaults.copy()
        param_dict[hyperparameter] = param

        # perform optimization
        start = time.time() #track the time
        results[param] = run_GA_custom_hyperparams(param_dict)
        comp_time = time.time()-start
        results[param] ={**results[param], 'time':comp_time}
    reversed_results = reverse_nested_dictionary(results)

    plot_hyperparams(reversed_results, parameter= hyperparameter,save=save)

def reverse_nested_dictionary(dict_to_reverse:dict, new_key:str='parameter')->dict:
    """
    must make sure all the sub dictionaries have the same layout
    :param dict_to_reverse: nested dictionary to reverse in the form of:
           `{key:{subkey:value, subkey2:value2, subkey3:value3}}`
    :param new_key: str: new key to save the old keys of the upperlevel dictionaries
    :return: reversed_dict: unnested dictionary where the keys of the subdictionary are the keys of upper dictionary
            `{new_key: list_of_keys, subkey:list_of_value, subkey2:list_of_value2, subkey3:list_of_value3}`
    """
    reversed_dict = {new_key: []}
    for key, sub_dict in dict_to_reverse.items():
        reversed_dict[new_key] += [key]
        for subkey,value in sub_dict.items():
            if subkey not in reversed_dict.keys():
                reversed_dict[subkey] = [value]
            else:
                reversed_dict[subkey] += [value]
    return reversed_dict


def plot_hyperparams(result_dict:dict,parameter:str = 'hyperparameter', save:bool=False):
    fig = plt.figure(layout="constrained")

    gs0 = GridSpec(1,2,figure=fig)
    gs_error_time = gs0[0,0]
    gs_kcats = gs0[0,1]

    gs = GridSpecFromSubplotSpec(2, 1,
                                height_ratios=[1, 1], hspace=0, wspace = 5,subplot_spec=gs_error_time)
    gs_error, gs_time = gs[0,0], gs[1,0]

    x_param = result_dict['parameter']

    ax_error = fig.add_subplot(gs_error)
    ax_error.plot(x_param,result_dict['error'])
    ax_error.set_title('Evolution of error')
    ax_error.set_xlabel(parameter)
    ax_error.set_ylabel('$R^{2}$')

    kcat1, kcat2,kcat3 = map(list, zip(*result_dict['kcats']))
    gs = GridSpecFromSubplotSpec(3, 1,
                                          height_ratios=[1, 1, 1], hspace=0, subplot_spec=gs_kcats)
    ax_kcat1 = fig.add_subplot(gs[0,0])
    ax_kcat2 = fig.add_subplot(gs[1,0])
    ax_kcat3 = fig.add_subplot(gs[2,0])

    ax_kcat1.plot(x_param, kcat1,color = 'darkred', label = 'E3')
    ax_kcat1.hlines(5, x_param[0], x_param[-1],color = 'darkred', linestyles='dashed')
    ax_kcat1.xaxis.set_visible(False)

    ax_kcat2.plot(x_param, kcat2, color = 'darkblue',label = 'E4')
    ax_kcat2.hlines(0.1, x_param[0], x_param[-1],color = 'darkblue', linestyles='dashed')
    ax_kcat2.xaxis.set_visible(False)

    ax_kcat3.plot(x_param, kcat3,color = 'purple', label = 'E5')
    ax_kcat3.hlines(0.25, x_param[0], x_param[-1], color = 'purple', linestyles='dashed')

    ax_kcat1.set_title('Evolution of kcat values')
    ax_kcat3.set_xlabel(parameter)
    ax_kcat2.set_ylabel('kcat value [$h^{-1}$]')

    ax_kcat1.set_ylim([0, 7])
    ax_kcat2.set_ylim([0, 2])
    ax_kcat3.set_ylim([0, 2])



    ax_time = fig.add_subplot(gs_time)
    ax_time.plot(x_param, result_dict['time'])
    ax_time.set_title('Computational time for optimization')
    ax_time.set_xlabel(parameter)
    ax_time.set_ylabel('time [s]')

    fig.set_figwidth(10)
    fig.set_figheight(7)
    os.chdir('../../Visualization/')
    if save: plt.savefig('Figures_hyperparams/'+parameter+'.png')
    else: plt.show()
    os.chdir('../Testing/Genetic_algorithm_tests/')


if __name__ == '__main__':
    #do simulations with different hyperparameters to see which one would be the best
    i=0
    while i<1:
        for param in hyperparameters.keys():
            hyperparam_range(param, hyperparameters[param], save =False)
            i+=1

