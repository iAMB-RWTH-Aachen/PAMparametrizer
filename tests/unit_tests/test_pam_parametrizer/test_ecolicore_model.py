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

def test_if_running_toy_model_in_pam_parametrizer_gives_correct_results():
    # Arrange
    sut = set_up_pamparametrizer(-10,-0.1)
    sut.validation_data._reactions_to_validate = RXNS_TO_VALIDATE
    sut._init_results_objects()
    bin_information = {1:[-10,-0.1, 1e-1]}
    reference_flux_data = run_pam_simulations(sut.pamodel)

    # Act
    sut.run_pamodel_simulations_in_bin(bin_information)
    flux_data = sut.parametrization_results.fluxes_df.drop('bin', axis =1)
    flux_data = flux_data.rename(columns = {'substrate': config.GLUCOSE_EXCHANGE_RXNID})
    flux_data = flux_data.sort_values(config.GLUCOSE_EXCHANGE_RXNID,ascending=True).reset_index(drop=True)

    # Assert
    assert pd.testing.assert_frame_equal(reference_flux_data, flux_data,
                                         check_dtype=False, check_exact=False, atol=1e-6) is None

def test_if_ecolicore_pam_old_params_has_better_rsquared_value_than_new_params():
    # Arrange
    pam_params_path = os.path.join('Data', 'proteinAllocationModel_iML1515_EnzymaticData_py.xls')
    old_params_pam = setup_ecolicore_pam(pam_data_file_path = pam_params_path)
    new_params_pam = setup_ecolicore_pam()


    # Act
    final_error_sut , sut, substrate_range = calculate_simulation_error_with_parametrizer_for_model(old_params_pam)
    sampled_data = adaptive_sampling(sut.validation_data.valid_data_df)


    # final_error_validation = evaluate_model_fitness_interpolation(toy_model=sut.pamodel,
    #                                                     substrate_rates=substrate_range,
    #                                                     validation_results = sut.validation_data.valid_data_df,
    #                                                     substrate_rxn = config.GLUCOSE_EXCHANGE_RXNID)
    final_error_validation = evaluate_model_fitness(model=sut.pamodel,
                                                        substrate_rates=sampled_data['EX_glc__D_e_ub'],#substrate_range,
                                                        validation_results = sampled_data,
                                                        substrate_rxn = config.GLUCOSE_EXCHANGE_RXNID)

    final_error_new_params, parametrizer_new, substrate_range = calculate_simulation_error_with_parametrizer_for_model(
        new_params_pam)

    # #plot results and save plot to see what is going on
    # fig, ax = sut.plot_valid_data()
    # parametrizer_new.plot_simulation(fig, ax)
    # sut.plot_simulation(fig, ax)
    # import matplotlib.pyplot as plt
    # plt.savefig('test.png')

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
    return result_df

def evaluate_model_fitness(model, validation_results:pd.DataFrame,
                               substrate_rates = [0.001, 0.091],
                               substrate_rxn:str = 'EX_glc__D_e') -> float:
    """
    Evaluate the fitness of the toymodel compared to the reference dataset generated using kcat_fwd = [1, 0.5, 5, 0.1, 0.25, 1.5]
    :return: float: error average difference of validation and result for the total of substrate uptake range and available reactiosn
    """
    simulation_results = run_pam_simulations(model, [abs(rate) for rate in substrate_rates])
    ref_data_rxn  = validation_results[validation_results[substrate_rxn+'_ub'].isin(substrate_rates)]
    error = []
    for rxn in RXNS_TO_VALIDATE:
        print(simulation_results[rxn])
        ref_data_rxn['simulation'] = simulation_results[rxn]
        # simulation mean
        data_average = ref_data_rxn[rxn].mean()
        # error: squared difference
        ref_data_rxn = ref_data_rxn.assign(error=lambda x: (x[rxn] - x['simulation']) ** 2)

        # calculate R^2:
        residual_ss = np.nansum(ref_data_rxn.error)
        total_ss = np.nansum([(data - data_average) ** 2 for data in ref_data_rxn[rxn]])
        # calculating r_squared is only feasible of the numerator and the denomenator are both nonzero
        if (residual_ss == 0) | (total_ss == 0):
            r_squared = 0
        else:
            r_squared = 1 - residual_ss / total_ss

        error += [r_squared]
    return np.nansum(error)/len(error)

# def sliding_window_r_squared(exp_x:list, exp_y:list, sim_x:list, sim_y:list, window_size:int) -> np.array:
#     """
#     Calculate the R^2 value between experimental and simulated data points within a sliding window.
#
#     Args:
#     - exp_x (array-like): X-coordinates of experimental data points.
#     - exp_y (array-like): Y-coordinates of experimental data points.
#     - sim_x (array-like): X-coordinates of simulated data points.
#     - sim_y (array-like): Y-coordinates of simulated data points.
#     - window_size (int): Size of the sliding window.
#
#     Returns:
#     - r_squared_values (array): squared error values computed for each sliding window.
#     """
#     r_squared_values = []
#     #get a pattern to find neighbouring x coordinated based on the sliding window size
#     pattern_list = list(range(-window_size // 2, window_size// 2 + 1))
#     for i in range(len(sim_x)):
#         # Find the indices of experimental data points within the sliding window
#         index = find_closest_datapoint(exp_x,abs(sim_x[i]))
#         window_indices = [index - diff for diff in pattern_list]
#         if len(window_indices) >= 2:
#             # Extract experimental and simulated data points within the sliding window
#             window_exp_x = exp_x[window_indices]
#             window_exp_y = exp_y[window_indices]
#             window_sim_y = np.interp(window_exp_x, sim_x, sim_y)
#             # Calculate R^2 value for the sliding window
#             slope, intercept, r_value, p_value, std_err = linregress(window_exp_y, window_sim_y)
#             simulation = [slope*x + intercept for x in window_exp_x]
#             # simulation mean
#             simulation_average = np.nanmean(simulation)
#             exp_average = np.mean(window_exp_y)
#             error = [(sim-abs(y))**2 for sim, y in zip(simulation, window_exp_y)]
#             # calculate R^2:
#             residual_ss = np.nansum(error)
#             total_ss = np.nansum([(data - exp_average) ** 2 for data in window_exp_y])
#
#             # calculating r_squared is only feasible of the numerator and the denomenator are both nonzero
#             if (residual_ss == 0) | (total_ss == 0):
#                 r_squared = 0
#             else:
#                 r_squared = 1 - residual_ss / total_ss
#             r_squared_values.append(r_squared)
#     return np.array(r_squared_values)

def interpolated_r_squared(exp_x:list, exp_y:list, sim_x:list, sim_y:list) -> np.array:
    """
    Calculate the R^2 value between experimental and simulated data points within a sliding window.

    Args:
    - exp_x (array-like): X-coordinates of experimental data points.
    - exp_y (array-like): Y-coordinates of experimental data points.
    - sim_x (array-like): X-coordinates of simulated data points.
    - sim_y (array-like): Y-coordinates of simulated data points.
    - window_size (int): Size of the sliding window.

    Returns:
    - r_squared_values (array): squared error values computed for each sliding window.
    """
    error = []
    for i in range(len(exp_x)):
        # Find the indices of experimental data points within the sliding window
        simulation = interpolate_datapoint(sim_x, sim_y, exp_x[i])
        error += [(simulation - abs(exp_y[i]))**2]
    # calculate R^2:
    exp_average = np.nanmean(exp_y)
    residual_ss = np.nansum(error)
    total_ss = np.nansum([(data - exp_average) ** 2 for data in exp_y])
    # calculating r_squared is only feasible of the numerator and the denomenator are both nonzero
    if (residual_ss == 0) | (total_ss == 0):
        r_squared = 0
    else:
        r_squared = 1 - residual_ss / total_ss

    return r_squared

def find_closest_datapoint(x_values, target_x):
    """
    Find the data point closest to a specific x value.

    Parameters:
    - x_values (array-like): Array of x-coordinates.
    - target_x (float): Target x value.

    Returns:
    - closest_index (int): Index of the data point closest to the target x value.
    """
    # Calculate absolute differences between target x value and each x-coordinate
    absolute_differences = np.abs(x_values - target_x)
    # Find index of the minimum absolute difference
    closest_index = np.argmin(absolute_differences)
    return closest_index

# def evaluate_model_fitness_sliding_window(toy_model, validation_results: pd.DataFrame,
#                            substrate_rates=[0.001, 0.091],
#                            substrate_rxn: str = 'EX_glc__D_e',
#                            window_size: int = 3) -> float:
#     """
#     Evaluate the fitness of the toy model compared to the reference dataset generated using kcat_fwd = [1, 0.5, 5, 0.1, 0.25, 1.5]
#     :return: float: error average difference of validation and result for the total of substrate uptake range and available reactions
#     """
#     simulation_results = run_pam_simulations(toy_model, substrate_rates)
#     error = []
#     for rxn in RXNS_TO_VALIDATE:
#         # Get experimental data points
#         exp_x = validation_results[substrate_rxn + '_ub']
#         exp_y = validation_results[rxn]
#         # Get simulated data points
#         sim_x = substrate_rates
#         sim_y = simulation_results[rxn].iloc[0:len(substrate_rates)]
#         # Calculate R^2 values using sliding window approach
#         r_squared_values = sliding_window_r_squared(exp_x, exp_y, sim_x, sim_y, window_size)
#         # Average R^2 values over sliding windows
#         avg_r_squared = np.mean(r_squared_values)
#         error.append(avg_r_squared)
#     return np.nanmean(error)

def evaluate_model_fitness_interpolation(toy_model,
                                         validation_results: pd.DataFrame,
                           substrate_rates=[0.001, 0.091],
                           substrate_rxn: str = 'EX_glc__D_e') -> float:
    """
    Evaluate the fitness of the toy model compared to the reference dataset generated using kcat_fwd = [1, 0.5, 5, 0.1, 0.25, 1.5]
    :return: float: error average difference of validation and result for the total of substrate uptake range and available reactions
    """
    simulation_results = run_pam_simulations(toy_model, substrate_rates)
    error = []
    for rxn in RXNS_TO_VALIDATE:
        # Get experimental data points
        exp_x = validation_results[substrate_rxn + '_ub']
        exp_y = validation_results[rxn]
        # Get simulated data points
        sim_x = substrate_rates
        sim_y = simulation_results[rxn].iloc[0:len(substrate_rates)]
        # Calculate R^2 values using sliding window approach
        r_squared = interpolated_r_squared(exp_x, exp_y, sim_x, sim_y)
        # Average R^2 values over sliding windows

        error.append(r_squared)
    return np.nanmean(error)

def interpolate_datapoint(x_values, y_values, x_target):
    """
    Interpolate a data point between two different x coordinates.

    Parameters:
    - x_values (array-like): Array of x-coordinates of the data points.
    - y_values (array-like): Array of y-coordinates of the data points.
    - x_target (float): Target x-coordinate to interpolate the data point.

    Returns:
    - interpolated_datapoint (float): Interpolated y-coordinate corresponding to the target x-coordinate.
    """
    # Find the indices of the two nearest x-coordinates
    idx_left = find_closest_datapoint(x_values, abs(x_target))
    idx_right = idx_left + 1
    if idx_right >= len(x_values)-1: idx_right = idx_left-1

    # Perform linear interpolation
    x_left, x_right = x_values[idx_left], x_values[idx_right]
    y_left, y_right = y_values[idx_left], y_values[idx_right]
    interpolated_datapoint = y_left + (y_right - y_left) * (abs(x_target) - x_left) / (x_right - x_left)
    # print('datapoint: ', x_target, 'closest simulated: ', x_left, x_right)
    # print('y simulated: ', y_left, y_right, 'y interpolated: ', interpolated_datapoint)

    return interpolated_datapoint

def calculate_simulation_error_with_parametrizer_for_model(pamodel):
    parametrizer = set_up_pamparametrizer(-10, -0.1)
    parametrizer.pamodel = pamodel
    parametrizer.validation_data._reactions_to_validate = RXNS_TO_VALIDATE
    parametrizer._init_results_objects()
    fluxes, substrate_range = parametrizer.run_simulations_to_plot()

    return parametrizer.calculate_final_error(fluxes, substrate_range), parametrizer, substrate_range