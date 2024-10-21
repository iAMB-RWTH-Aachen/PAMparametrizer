from itertools import combinations
import pandas as pd
import os
from datetime import date

from Scripts.i2_parametrization.pam_parametrizer_performance_analysis import (get_statistics_from_df,
                                                                   save_pam_parametrizer_results_to_df)
from Modules.utils.error_calculation import nanaverage


from Scripts.i2_parametrization.pam_parametrizer_iML1515 import set_up_pamparametrizer, set_up_validation_data

RESULT_FILE_PATH = os.path.join('Results', f'pam_parametrizer_statistics_leave_out_csource_{date.today()}.xlsx')


def run_parametrization_workflow(iteration, iterations,
                                 processes, leave_out_csource,
                                 gene_flow_events, num_kcats_to_mutate,
                                 best_individual_df, computational_performance_df,
                                 rsquared_df = None,
                                 min_substrate_uptake = -11, max_substrate_uptake = -0.1):
    print('\n\n###################################################################################')
    print('starting with iteration number ', iteration, ' out of ', iterations, ' iterations with data for',
          " and ".join(leave_out_csource))
    print('------------------------------------------------------------------------------------------------')
    all_csources = ['Glycerol', 'Glucose', 'Acetate', 'Pyruvate', 'Gluconate', 'Succinate', 'Galactose', 'Fructose']
    # csources =[c for c in all_csources if c not in leave_out_csource]
    csources = ['Glucose']+leave_out_csource

    parametrizer = set_up_pamparametrizer(min_substrate_uptake, max_substrate_uptake, processes=processes,
                                              gene_flow_events=gene_flow_events,
                                              filename_extension= f"_{'_'.join(leave_out_csource)}",
                                              num_kcats_to_mutate = num_kcats_to_mutate,
                                          kcat_increase_factor = 3,
                                          c_sources = csources)

    #need to reset best individual and computational performance df
    parametrizer.parametrization_results.best_individuals = pd.DataFrame(columns=['run_id', 'enzyme_id',
                                                                                  'direction', 'rxn_id',
                                                                                  'kcat[s-1]', 'ga_error'])
    parametrizer.parametrization_results.computational_time = pd.DataFrame(columns=['run_id', 'time_s', 'time_h'])

    parametrizer.run(binned='False')

    results_file_path = os.path.join('Results', f"pam_parametrizer_diagnostics_{'_'.join(leave_out_csource)}.xlsx")
    best_individual_df, computational_performance_df = save_pam_parametrizer_results_to_df(iteration,'_'.join(leave_out_csource),
                                                   best_individual_df,computational_performance_df,
                                                   results_file_path)
    rsquared_df = get_r_squared_for_all_csources(parametrizer, iteration, {'_'.join(leave_out_csource)}, all_csources,
                                                 rsquared_df)

    return best_individual_df, computational_performance_df, rsquared_df

def analyse_parametrizer_performance():
    min_substrate = -11
    max_substrate = -0.1
    iterations = 4
    processes = 2
    gene_flow_events = processes
    all_csources = ['Glycerol', 'Glucose', 'Acetate', 'Pyruvate', 'Gluconate', 'Succinate', 'Galactose', 'Fructose']
    csource_combinations= get_csource_combinations(all_csources)

    # 0. Initialize result dataframes
    best_individual_df = pd.DataFrame(columns=['iteration', 'binned',
                                               'run_id', 'enzyme_id', 'rxn_id', 'kcat[s-1]',
                                               'ga_error', 'r_squared'])
    computational_performance_df = pd.DataFrame(columns=['iteration', 'binned', 'run_id', 'time_s', 'time_h'])
    rsquared_df = None

    # 1. Run different configuations of the toy model parametrization for n iterations and save the results
    for csources_leavout in csource_combinations:
        for iteration in range(iterations):
            best_individual_df, computational_performance_df, rsquared_df = run_parametrization_workflow(iteration+1, iterations,
                                     processes,csources_leavout,
                                     gene_flow_events, 5,
                                     best_individual_df, computational_performance_df, rsquared_df,
                                     min_substrate_uptake=min_substrate, max_substrate_uptake=max_substrate)

    # 2. Calculate mean error and stdev for each method and for a single iteration
    error_per_iteration_config = get_statistics_from_df(best_individual_df,
                                                        group_by=['iteration', 'binned'],
                                                        columns=['ga_error', 'r_squared'])
    error_per_config = get_statistics_from_df(best_individual_df,
                                              group_by=['binned'],
                                              columns=['ga_error', 'r_squared'])
    computational_performance_per_config = get_statistics_from_df(computational_performance_df,
                                                                  group_by=['binned'],
                                                                  columns=['time_s'])
    # 3. save the results to excel
    with pd.ExcelWriter(RESULT_FILE_PATH) as writer:
        # Write each DataFrame to a specific sheet
        best_individual_df.to_excel(writer, sheet_name='Best_Individuals', index=False)
        computational_performance_df.to_excel(writer, sheet_name='Computational_Time',
                                              index=False)
        rsquared_df.to_excel(writer, sheet_name='R2_per_csource', index=False)
        error_per_config.to_excel(writer, sheet_name='stats_error_per_config', index=False)
        error_per_iteration_config.to_excel(writer, sheet_name='stats_error_per_iteration', index=False)
        computational_performance_per_config.to_excel(writer, sheet_name='stats_computation_time', index=False)

def get_csource_combinations(csources):
    # Generate all combinations of two strings
    combinations_of_two = list(combinations(csources, 2))

    # Convert the tuples to lists
    combinations_of_two = [list(combo) for combo in combinations_of_two]

    # Combine individual strings with their combinations
    result = [[item] for item in csources] + combinations_of_two

    return result

def get_r_squared_for_all_csources(parametrizer, iteration:int, configuration:str, csources:list[str],
                                   rsquared_df:pd.DataFrame = None):
    if rsquared_df is None:
        rsquared_df = pd.DataFrame(columns=["iteration", "configuration"]+csources)

    parametrizer.validation_data = set_up_validation_data(csources)
    run_simulations_to_validate(parametrizer)

    error = [iteration, configuration]
    for substrate_uptake_id in csources:
        valid_data = parametrizer.validation_data.get_by_id(substrate_uptake_id)
        reactions_to_validate = valid_data._reactions_to_validate
        validation_df = valid_data.sampled_valid_data  # .apply(lambda x: x.abs() if x.dtype!='object' else x)

        # self._init_validation_df(bin_information=[self.min_substrate_uptake_rate, self.max_substrate_uptake_rate])
        error += [nanaverage(parametrizer._calculate_error_for_reactions(substrate_uptake_id=substrate_uptake_id,
                                                                 validation_df=validation_df,
                                                                 reactions_to_validate=reactions_to_validate,
                                                                 bin_id="final"))]
    rsquared_df.loc[len(rsquared_df)] = error
    return rsquared_df

def run_simulations_to_validate(parametrizer):
    for valid_data in parametrizer.validation_data:
        sampled_data = valid_data.sampled_valid_data
        substrate_uptake_id = valid_data.id
        substrate_rates_for_error = sampled_data[substrate_uptake_id + "_ub"].to_list()

        fluxes, substrate_rates = parametrizer.run_simulations_to_plot(substrate_uptake_id,
                                                               substrate_rates=substrate_rates_for_error,
                                                                       sensitivity=False)

        for simulation_result, substrate_rate in zip(fluxes, substrate_rates):
            parametrizer.parametrization_results.add_fluxes_from_fluxdict(flux_dict=simulation_result,
                                                                  bin_id="final",
                                                                  substrate_reaction_id=substrate_uptake_id,
                                                                  substrate_uptake_rate=substrate_rate,
                                                                  fluxes_abs=False)


if __name__ == '__main__':
    analyse_parametrizer_performance()