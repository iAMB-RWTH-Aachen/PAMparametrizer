import numpy as np
import time
from multiprocessing import Pool
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from typing import Iterable


#ignoring infeasible warnings for cleaner output
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

os.chdir('..')
from Scripts.genetic_algorithm_analysis import get_kcat_error_from_ga, reverse_nested_dictionary

from Scripts.Testing.Genetic_algorithm_tests.toy_model import init_toy_parametrization_ga


sys.path.append('../../i3_analysis')

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
    return get_kcat_error_from_ga(ga)

def hyperparam_range(hyperparameter: str, range: Iterable[float], save:bool=False):
    final_results = []
    i=0
    while i<=20:
        results = {}
        i+= 1
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
        final_results += [reversed_results]

    plot_hyperparams(final_results, parameter= hyperparameter,save=save)

def plot_hyperparams(result_list:list, parameter:str = 'hyperparameter', save:bool=False):
    fig = plt.figure(layout="constrained")

    gs0 = GridSpec(1,2,figure=fig)
    gs_error_time = gs0[0,0]
    gs_kcats = gs0[0,1]

    gs = GridSpecFromSubplotSpec(2, 1,
                                height_ratios=[1, 1], hspace=0, wspace = 5,subplot_spec=gs_error_time)
    gs_error, gs_time = gs[0,0], gs[1,0]
    gs = GridSpecFromSubplotSpec(3, 1,
                                 height_ratios=[1, 1, 1], hspace=0, subplot_spec=gs_kcats)
    ax_kcat1 = fig.add_subplot(gs[0, 0])
    ax_kcat2 = fig.add_subplot(gs[1, 0])
    ax_kcat3 = fig.add_subplot(gs[2, 0])
    x_param = result_list[0]['parameter']
    ax_error = fig.add_subplot(gs_error)
    ax_time = fig.add_subplot(gs_time)

    for result_dict in result_list:
        ax_error.plot(x_param, result_dict['error'])
        ax_error.set_title('Evolution of error')
        ax_error.set_xlabel(parameter)
        ax_error.set_ylabel('$R^{2}$')

        kcat1, kcat2,kcat3 = map(list, zip(*result_dict['kcats']))

        ax_kcat1.plot(x_param, kcat1,color = 'darkred', label = 'E3')
        ax_kcat1.hlines(5, x_param[0], x_param[-1],color = 'darkred', linestyles='dashed')
        ax_kcat1.xaxis.set_visible(False)

        ax_kcat2.plot(x_param, kcat2, color = 'darkblue',label = 'E4')
        ax_kcat2.hlines(0.1, x_param[0], x_param[-1],color = 'darkblue', linestyles='dashed')
        ax_kcat2.xaxis.set_visible(False)

        ax_kcat3.plot(x_param, kcat3,color = 'purple', label = 'E5')
        ax_kcat3.hlines(0.25, x_param[0], x_param[-1], color = 'purple', linestyles='dashed')

        ax_time.plot(x_param, result_dict['time'])

    ax_kcat1.set_title('Evolution of kcat values')
    ax_kcat3.set_xlabel(parameter)
    ax_kcat2.set_ylabel('kcat value [$h^{-1}$]')

    ax_kcat1.set_ylim([0, 7])
    ax_kcat2.set_ylim([0, 2])
    ax_kcat3.set_ylim([0, 2])

    ax_time.set_title('Computational time for optimization')
    ax_time.set_xlabel(parameter)
    ax_time.set_ylabel('time [s]')

    fig.set_figwidth(10)
    fig.set_figheight(7)
    if save: plt.savefig('Scripts/i3_analysis/Figures_hyperparams/'+parameter+'.png')
    else: plt.show()



if __name__ == '__main__':
    #do simulations with different hyperparameters to see which one would be the best

    for param in hyperparameters.keys():
        hyperparam_range(param, hyperparameters[param], save =True)
        break


