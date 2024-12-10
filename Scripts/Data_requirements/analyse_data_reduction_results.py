from typing import Iterable
import numpy as np
import pandas as pd
import os

from Scripts.i3_analysis.PAMparametrizer_compare_alternative_solutions import set_up_ecoli_pam_parametrizer_and_get_substrate_uptake_rates

from Modules.utils.pam_generation import create_pamodel_from_diagnostics_file

def run_simulations_and_calculate_error(parametrizer, substrate_rates) -> float:
    fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                     substrate_rates=substrate_rates,
                                                     sensitivity=False)
    fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                     substrate_rates=substrate_rates,
                                                     sensitivity=False)

    for flux, rate in zip(fluxes, substrate_rates):
        parametrizer.parametrization_results.add_fluxes_from_fluxdict(flux_dict=flux,
                                                                      bin_id='final',
                                                                      substrate_reaction_id=parametrizer.substrate_uptake_id,
                                                                      substrate_uptake_rate=rate,
                                                                      fluxes_abs=False)
    error = parametrizer.calculate_final_error()
    return error

def get_error_for_parametrization_experiment(parametrizer, best_indiv_file_path:str, substrate_rates: Iterable) -> float:
    pamodel = parametrizer.pamodel.copy(copy_with_pickle=True)
    parametrizer.pamodel = create_pamodel_from_diagnostics_file(best_indiv_file_path, pamodel)
    return run_simulations_and_calculate_error(parametrizer, substrate_rates)

def compute_final_error_on_full_dataset_for_all_experiments(base_file_path: str,
                                                            datasizes: Iterable,
                                                            num_replicates:float) -> pd.DataFrame:
    parametrizer, substrate_rates = set_up_ecoli_pam_parametrizer_and_get_substrate_uptake_rates()
    #ensure all the errors are calculated based pn the same dataset
    parametrizer.validation_data.get_by_id('EX_glc__D_e').sampled_valid_data = parametrizer.validation_data.get_by_id('EX_glc__D_e').valid_data

    final_errors = pd.DataFrame(columns = ['perc_data', 'sample', 'final_error'])

    for datasize in datasizes:
        print('\n------------------------------------------------------')
        print(f"Analyzing parametrization with {datasize}% of the total amount of data to train")
        for sample in range(1,num_replicates+1):
            print(f'\tReplicate {sample}')
            file_path = f'{base_file_path}{datasize}_{sample}.xlsx'
            error = get_error_for_parametrization_experiment(parametrizer, file_path, substrate_rates)
            final_errors.loc[len(final_errors)] = [datasize, sample, error]

if __name__ == '__main__':
    diagnostic_file_path_base = os.path.join('Results', 'data_reduction_results', 'diagnostics', 'pam_parametrizer_diagnostics_datareduc_')
    final_errors = compute_final_error_on_full_dataset_for_all_experiments(diagnostic_file_path_base,
                                                                           datasizes= np.arange(10,80,10),
                                                                           num_replicates=3)
    print(final_errors)

