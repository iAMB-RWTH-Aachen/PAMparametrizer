from Scripts.Testing.pam_parametrizer_ecolicore import set_up_pamparametrizer
from Scripts.pam_generation import setup_ecolicore_pam

from scipy.stats import linregress
import pytest
from PAModelpy.configuration import Config
import pandas as pd
import numpy as np
import os

config = Config()
config.reset()
RXNS_TO_VALIDATE = [config.ACETATE_EXCRETION_RXNID, config.CO2_EXHANGE_RXNID, config.OXYGEN_UPTAKE_RXNID, 'BIOMASS_Ecoli_core_w_GAM']

# def test_if_running_toy_model_in_pam_parametrizer_gives_correct_results():
#     # Arrange
#     sut = set_up_pamparametrizer(-10,-0.1)
#     sut.validation_data._reactions_to_validate = RXNS_TO_VALIDATE
#     sut._init_results_objects()
#     bin_information = {1:[-10,-0.1, 1e-1]}
#     substrate_range = [abs(rate) for rate in set(sut.validation_data.valid_data_df[sut.substrate_uptake_id+'_ub'].to_list())]
#     print(substrate_range)
#     reference_flux_data = run_pam_simulations(sut.pamodel, substrate_rates=substrate_range)
#
#     # Act
#     sut.run_pamodel_simulations_in_bin(bin_information)
#     flux_data = sut.parametrization_results.fluxes_df.drop('bin', axis =1)
#     flux_data = flux_data.rename(columns = {'substrate': config.GLUCOSE_EXCHANGE_RXNID}).drop_duplicates()
#     flux_data = flux_data.sort_values(config.GLUCOSE_EXCHANGE_RXNID,ascending=True).reset_index(drop=True)
#     for i in range(len(reference_flux_data)):
#         if not reference_flux_data.iloc[i].equals(flux_data.iloc[i]):
#             print(reference_flux_data.iloc[i], flux_data.iloc[i])
#     print(reference_flux_data, flux_data)
#     # Assert
#     assert pd.testing.assert_frame_equal(reference_flux_data, flux_data,
#                                          check_dtype=False, check_exact=False, atol=1e-6) is None

def test_if_ecolicore_pam_old_params_has_better_rsquared_value_than_new_params():
    # Arrange
    pam_params_path = os.path.join('Data', 'proteinAllocationModel_iML1515_EnzymaticData_py.xls')
    old_params_pam = setup_ecolicore_pam(pam_data_file_path = pam_params_path)
    new_params_pam = setup_ecolicore_pam()


    # Act
    final_error_sut , sut, substrate_range = calculate_simulation_error_with_parametrizer_for_model(old_params_pam)

    final_error_new_params, parametrizer_new, substrate_range = calculate_simulation_error_with_parametrizer_for_model(
        new_params_pam)
    sampled_data = sut.validation_data.get_by_id('EX_glc__D_e').sampled_valid_data

    final_error_validation = evaluate_model_fitness(model=sut.pamodel,
                                                    substrate_rates=sampled_data['EX_glc__D_e_ub'],#substrate_range,
                                                    validation_results = sampled_data,
                                                    substrate_rxn = config.GLUCOSE_EXCHANGE_RXNID)

    # assert
    # r^2 for new params should be smaller than r^2 for parametrized model
    assert final_error_new_params < final_error_sut
    assert final_error_validation == pytest.approx(final_error_sut, 1e-2)


########################################################################################################################
# HELPER FUNCTION
########################################################################################################################
def adaptive_sampling(exp_data:pd.DataFrame,substrate_rxn: str = 'EX_glc__D_e_ub', num_samples=10, min_density=5, max_density=20):
    """
    Perform adaptive sampling based on data point distribution and variability.

    Parameters:
    - x_exp (array-like): X-coordinates of experimental data.
    - y_exp (array-like): Y-coordinates of experimental data.
    - num_samples (int): Number of samples to generate.
    - min_density (int): Minimum density of samples per unit distance.
    - max_density (int): Maximum density of samples per unit distance.

    Returns:
    - sampled_x (array): Sampled x-coordinates.
    """
    # Calculate distances between consecutive data points
    distances = np.diff(exp_data[substrate_rxn])

    # Compute local data point density as reciprocal of distance
    densities = 1 / distances

    # Normalize densities to [min_density, max_density]
    densities = np.clip(densities, min_density, max_density)

    # Calculate sampling probabilities based on normalized densities
    probabilities = densities / np.sum(densities)

    # Sample indices with replacement based on probabilities
    sampled_indices = np.random.choice(np.arange(len(exp_data) - 1), size=num_samples, p=probabilities)

    # Sample x-coordinates from selected indices
    sampled = exp_data.loc[sampled_indices] #+ np.random.uniform(low=0, high=distances[sampled_indices])

    return sampled


def run_pam_simulations(pamodel, substrate_rates = np.arange(0.1,10.1,0.1)):
    result_df = pd.DataFrame(columns=[config.GLUCOSE_EXCHANGE_RXNID]+RXNS_TO_VALIDATE)

    for substrate in substrate_rates:
        pamodel.change_reaction_bounds(rxn_id=config.GLUCOSE_EXCHANGE_RXNID,
                                       lower_bound=-substrate, upper_bound=0)
        print('Running simulations with ', substrate, 'mmol/g_cdw/h of substrate going into the system')
        pamodel.optimize()
        if pamodel.solver.status == 'optimal' and pamodel.objective.value > 0:
            results_row = []
            for rxn in RXNS_TO_VALIDATE:
                results_row += [abs(pamodel.reactions.get_by_id(rxn).flux)]
            result_df.loc[len(result_df)] = [substrate] + results_row
    result_df = result_df.sort_values(by =config.GLUCOSE_EXCHANGE_RXNID, ascending = True).reset_index(drop=True)
    return result_df.drop_duplicates()

def evaluate_model_fitness(model, validation_results:pd.DataFrame,
                               substrate_rates = [0.001, 0.091],
                               substrate_rxn:str = 'EX_glc__D_e') -> float:
    """
    Evaluate the fitness of the toymodel compared to the reference dataset generated using kcat_fwd = [1, 0.5, 5, 0.1, 0.25, 1.5]
    :return: float: error average difference of validation and result for the total of substrate uptake range and available reactiosn
    """
    simulation_results = run_pam_simulations(model, [abs(rate) for rate in substrate_rates])
    validation_df  = validation_results[validation_results[substrate_rxn+'_ub'].isin(substrate_rates)].apply(lambda x: x.abs() if x.dtype!='object' else x)
    # ref_data_rxn.update(ref_data_rxn.select_dtypes(include=[np.number]).abs())
    error = []
    for rxn in RXNS_TO_VALIDATE:
        validation_data = validation_df.dropna(axis=0, subset=rxn)
        # if there are no reference data points, continue to the next reaction
        if len(validation_data) == 0:
            continue
        simulated_data = pd.DataFrame({substrate_rxn+'_ub': [abs(flux) for flux in simulation_results[substrate_rxn]],
                                       'simulation': simulation_results[rxn] })
        ref_data_rxn = pd.merge(validation_data,simulated_data, on=substrate_rxn+'_ub', how='inner')
        # error: squared difference
        ref_data_rxn = ref_data_rxn.assign(error=lambda x: (x[rxn] - x['simulation']) ** 2)
        # calculate R^2:
        data_average = np.nanmean(ref_data_rxn[rxn])
        residual_ss = np.nansum(ref_data_rxn.error)
        total_ss = np.nansum([(data - data_average) ** 2 for data in ref_data_rxn[rxn]])
        # calculating r_squared is only feasible of the numerator and the denomenator are both nonzero
        if (residual_ss == 0) | (total_ss == 0):
            r_squared = 0
        else:
            r_squared = 1 - residual_ss / total_ss

        error += [r_squared]
    return np.nanmean(error)

def calculate_simulation_error_with_parametrizer_for_model(pamodel):
    parametrizer = set_up_pamparametrizer(-10, -0.1)
    parametrizer.pamodel = pamodel
    parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_validate = RXNS_TO_VALIDATE
    parametrizer._init_results_objects()
    parametrizer._init_validation_df(bin_information=[-10, -0.1])
    substrate_rates = parametrizer.validation_data.get_by_id('EX_glc__D_e').sampled_valid_data['EX_glc__D_e_ub']
    fluxes, substrate_range = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                                   substrate_rates=substrate_rates)

    #save flux data
    # check if there is input data (only used in unit testing)
    for simulation_result, substrate_rate in zip(fluxes, substrate_rates):
        parametrizer.parametrization_results.add_fluxes_from_fluxdict(flux_dict=simulation_result,
                                                                      bin_id='final',
                                                                      substrate_reaction_id='EX_glc__D_e',
                                                                      substrate_uptake_rate=substrate_rate,
                                                                      fluxes_abs=False)

    return parametrizer.calculate_final_error(), parametrizer, substrate_range