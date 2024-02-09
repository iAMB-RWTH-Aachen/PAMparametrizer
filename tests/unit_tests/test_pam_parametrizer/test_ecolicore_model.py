from Scripts.Testing.pam_parametrizer_ecolicore import set_up_pamparametrizer
from PAModelpy.configuration import Config
import pandas as pd
import numpy as np

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