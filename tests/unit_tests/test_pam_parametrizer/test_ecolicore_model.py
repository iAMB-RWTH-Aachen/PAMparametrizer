from Scripts.i2_parametrization.pam_parametrizer_ecolicore import set_up_pamparametrizer
from Scripts.i2_parametrization.pam_parametrizer_iML1515 import set_up_pamparametrizer as set_up_pamparametrizer_iml1515
from Scripts.pam_generation_uniprot_id import setup_ecolicore_pam
from Scripts.pam_generation import setup_ecolicore_pam as setup_ecolicore_pam_ec


from tests.unit_tests.test_pam_parametrizer.test_pam_parametrizer import (save_simulated_fluxes_in_pamparametrizer_for_different_carbon_sources,
                                                                          get_kcat_values_from_model)

import pytest
from PAModelpy.configuration import Config
import pandas as pd
import numpy as np
import os

config = Config()
config.reset()
RXNS_TO_VALIDATE = [config.ACETATE_EXCRETION_RXNID, config.CO2_EXHANGE_RXNID, config.OXYGEN_UPTAKE_RXNID, 'BIOMASS_Ecoli_core_w_GAM']

#IMPORTANT:TESTS ARE DEPRICATED. NEED TO BE UPDATED ASAP TODO
# def test_if_ecolicore_pam_old_params_has_better_rsquared_value_than_gem():
#     # Arrange
#     pam_params_path = os.path.join('Data', 'proteinAllocationModel_iML1515_EnzymaticData_py.xls')
#
#     old_params_pam = setup_ecolicore_pam_ec(pam_data_file_path = pam_params_path)
#     gem = setup_ecolicore_pam_ec(pam_data_file_path = pam_params_path, total_protein=1e3) # large protein pool, thus acts as gem
#
#
#     # Act
#     final_error_sut ,sut, substrate_range = calculate_simulation_error_with_parametrizer_for_model(old_params_pam)
#
#     final_error_gem, parametrizer_gem, substrate_range = calculate_simulation_error_with_parametrizer_for_model(
#         gem)
#     sampled_data = sut.validation_data.get_by_id('EX_glc__D_e').sampled_valid_data
#
#     final_error_validation = evaluate_model_fitness(model=sut._pamodel,
#                                                     substrate_rates=sampled_data['EX_glc__D_e_ub'],  #substrate_range,
#                                                     validation_results = sampled_data,
#                                                     substrate_rxn = config.GLUCOSE_EXCHANGE_RXNID)
#
#     # assert
#     # r^2 for gem should be smaller than r^2 for parametrized model
#     assert final_error_gem < final_error_sut
#     assert final_error_validation == pytest.approx(final_error_sut, 1e-2)
#
# def test_if_simulation_error_for_multiple_carbon_sources_of_parametrizer_is_same_as_for_genetic_algorithm():
#     #setup PAMparametrizer
#     sut_param = set_up_pamparametrizer(-11, -0.1,other_csources = True)
#     sut_param.calculate_translational_sector_for_multiple_csources()
#     transl_sector_config = {vd.id: vd.translational_sector_config for vd in sut_param.validation_data}
#
#     sut_param._init_validation_df(substrate_uptake_ids = sut_param.substrate_uptake_ids)
#     sut_param._init_results_objects()
#
#     sut_param = save_simulated_fluxes_in_pamparametrizer_for_different_carbon_sources(sut_param)
#     #get substrate uptake rates to calculate the error for
#     substrate_rates = {vd.id: vd.sampled_valid_data[vd.id+'_ub'].values for vd in sut_param.validation_data}
#
#     # set up genetic algorithm
#     sut_ga = sut_param._init_genetic_algorithm(substrate_uptake_rates = substrate_rates,
#                                                enzymes_to_evaluate = {
#                                                    'P0A9Q7': {
#                                                        'reaction':'ALCD2x', 'kcats': {'b': 1/(1.6088* 3600*1e-6),
#                                                                                       'f': 1/(11.6366* 3600*1e-6)},
#                                                        'sensitivity': 0.1},
#                                                    'P08200': {
#                                                        'reaction': 'ICDHyr', 'kcats': {'f': 1/(10.6819* 3600*1e-6),
#                                                                                        'b': 1/(10.6819* 3600*1e-6)},
#                                                        'sensitivity': 0.1},
#                                                    'P0ABH7': {
#                                                        'reaction':'CS', 'kcats': {'f':1/(101.8665* 3600*1e-6)},
#                                                        'sensitivity': 0.1}
#                                                },
#                                                translational_sector_config=transl_sector_config,
#                                                filename_extension = '')
#     toolbox = sut_ga._init_deap_toolbox()
#     toy_ga = sut_ga.ga
#     population = toolbox.population(n=3)
#     population[0].kcat_list = [1/(1.6088* 3600*1e-6),1/(11.6366* 3600*1e-6),1/(10.6819* 3600*1e-6),
#                                1/(10.6819* 3600*1e-6),1/(101.8665* 3600*1e-6)]
#
#     # Act
#     #pam parametrizer
#     r_squared_param = sut_param.calculate_final_error()
#     #genetic algorithm
#     population = toy_ga.evaluate_pop(population, toolbox)
#     individual_to_evaluate = population[0]
#     r_squared_ga = individual_to_evaluate.fitness._wsum()
#
#     # Assert
#     assert r_squared_ga == pytest.approx(r_squared_param, abs = 1e-3)
#
# #
# def test_pam_parametrizer_configures_translational_sector_correctly():
#     # Arrange
#     sut = set_up_pamparametrizer_iml1515(-11, -0.1)
#     #reset translational sector
#     sut.validation_data.get_by_id('EX_glc__D_e').translational_sector_config = None
#     #get the reference parameters
#     tps_0 = sut._pamodel.sectors.get_by_id('TranslationalProteinSector').tps_0[0]
#     tps_mu = sut._pamodel.sectors.get_by_id('TranslationalProteinSector').tps_mu[0]
#
#     # Apply
#     sut.calculate_translational_sector_for_multiple_csources()
#     transl_sector_config_sut = sut.validation_data.get_by_id('EX_glc__D_e').translational_sector_config
#
#     # Assert
#     assert transl_sector_config_sut['intercept'] == pytest.approx(tps_0, rel=1e-1)
#     assert transl_sector_config_sut['slope'] == pytest.approx(tps_mu, rel=1e-1)
#
# def test_pam_parametrizer_changes_kcats_same_way_as_genetic_algorithm():
#     # Arrange
#     kcat = 5
#     enzyme_id = 'P26616'
#     reaction_id = 'ME1'
#     kcat_dict = {reaction_id: {'f': kcat}}
#     enzymes_to_evaluate = {enzyme_id: {
#                         'reaction': reaction_id,
#                         'kcats': {'f':kcat},
#                         'sensitivity': 0.1}}
#     substrate_rates = {'EX_glc_D__e':[-6,-5]}
#
#     #set up parametrizer and genetic algorithm objects
#     sut = set_up_pamparametrizer(-11, -0.1)
#     ga = sut._init_genetic_algorithm(substrate_rates,
#                                      enzymes_to_evaluate,
#                                      translational_sector_config= {'EX_glc__D_e':
#                                                                        sut.validation_data.EX_glc__D_e.translational_sector_config},
#                                      filename_extension=''
#                                      )
#     #for the genetic algorithm we need a dummy individual to change the kcat
#     toolbox = ga._init_deap_toolbox()
#     population = ga.ga.init_pop(toolbox, ga.population_size, True)
#     individual = population[0]
#     individual.kcat_list = [(kcat*3600*1e-6)]
#
#     # Act
#     sut._change_kcat_value_for_enzyme(enzyme_id=enzyme_id, kcat_dict=kcat_dict)
#     kcat_model_sut = get_kcat_values_from_model(sut._pamodel, enzyme_ids= [enzyme_id], reaction_names= [f"CE_{reaction_id}_{enzyme_id}"])[0]
#
#     ga.FitEval._change_kcat_values_for_individual(individual)
#     kcat_model_ga = get_kcat_values_from_model(ga.FitEval.model, enzyme_ids= [enzyme_id], reaction_names= [f"CE_{reaction_id}_{enzyme_id}"])[0]
#
#     # Assert
#     assert kcat_model_sut == pytest.approx(kcat_model_ga, abs = 1e-3)
#     assert kcat == pytest.approx(kcat_model_sut, abs = 1e-3)
#     assert kcat == pytest.approx(kcat_model_ga, abs = 1e-3)


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

def calculate_simulation_error_with_parametrizer_for_model(pamodel, other_csources = False):
    parametrizer = set_up_pamparametrizer(-10, -0.1, other_csources)
    parametrizer._pamodel = pamodel
    parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_validate = RXNS_TO_VALIDATE
    parametrizer._init_results_objects()
    parametrizer._init_validation_df(bin_information=[-10, -0.1])
    substrate_rates = parametrizer.validation_data.get_by_id('EX_glc__D_e').sampled_valid_data['EX_glc__D_e_ub']
    fluxes, substrate_range = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                                   substrate_rates=substrate_rates)



    #save flux data
    # check if there is input data (only used in unit testing)
    for simulation_result, substrate_rate in zip(fluxes, substrate_range.values()):
        parametrizer.parametrization_results.add_fluxes_from_fluxdict(flux_dict=simulation_result,
                                                                      bin_id='final',
                                                                      substrate_reaction_id='EX_glc__D_e',
                                                                      substrate_uptake_rate=substrate_rate[0],
                                                                      fluxes_abs=False)

    return parametrizer.calculate_final_error(), parametrizer, substrate_range.keys()

