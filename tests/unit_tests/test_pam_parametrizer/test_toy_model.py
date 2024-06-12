import pandas as pd

from Scripts.Testing.pam_parametrizer_toy_model import set_up_pamparametrizer

FINAL_ENZYMES2KCAT = {'E1':{'R1':{'f': 1/(3600*1e-6), 'b':1/(3600*1e-6)}},
                          'E2': {'R2': {'f': 0.5 / (3600 * 1e-6), 'b': 0.5 / (3600 * 1e-6)}},
                          'E3':{'R3':{'f': 5/(3600*1e-6), 'b':5/(3600*1e-6)}},
                      'E4':{'R4':{'f': 0.1/(3600*1e-6), 'b':0.1/(3600*1e-6)}},
                      'E5':{'R5':{'f': 0.25/(3600*1e-6), 'b':0.25/(3600*1e-6)}},
                          'E6':{'R6':{'f': 1.5/(3600*1e-6), 'b':1.5/(3600*1e-6)}}}

def test_if_toy_model_parameters_in_pam_parametrizer_are_set_correctly():
    # Arrange
    sut = set_up_pamparametrizer(0.001, 0.1, kcat_fwd = [1, 0.5, 5, 0.1, 0.25, 1.5])
    # change_kcat_to_expected_outcome(sut)

    # Assert
    for enzyme, rxn2kcat_dict in FINAL_ENZYMES2KCAT.items():
        kcat_to_validate = sut.pamodel.enzymes.get_by_id(enzyme).rxn2kcat
        assert all(v == kcat_to_validate[k] for k,v in rxn2kcat_dict.items()) and len(rxn2kcat_dict)==len(kcat_to_validate)

def test_if_running_toy_model_in_pam_parametrizer_gives_correct_results():
    # Arrange
    sut = set_up_pamparametrizer(0.001, 0.1, kcat_fwd = [1, 0.5, 5, 0.1, 0.25, 1.5])
    sut.validation_data.get_by_id('R1')._reactions_to_validate = ['R1', 'R7', 'R8', 'R9']
    sut._init_results_objects()
    bin_id = 1
    bin_information = [0.001,0.1, 1e-2]
    # change_kcat_to_expected_outcome(sut)
    reference_flux_data = sut.validation_data.get_by_id('R1').valid_data.iloc[:,1:]

    # Act
    sut.run_pamodel_simulations_in_bin(bin_id = bin_id, bin_information=bin_information, substrate_uptake_reaction='R1')
    flux_data = sut.parametrization_results.flux_results.get_by_id('R1').fluxes_df.drop('bin', axis =1)
    flux_data = flux_data.rename(columns = {'substrate': 'R1_ub'})

    # Assert
    assert pd.testing.assert_frame_equal(reference_flux_data, flux_data,
                                         check_dtype=False,check_exact=False, atol=1e-3) is None


##########################################################################################################################
# HELPER FUNCTIONS
##########################################################################################################################
def change_kcat_to_expected_outcome(pam_parametrizer):
    for enz_id, kcat_dict in FINAL_ENZYMES2KCAT.items():
        pam_parametrizer.pamodel.change_kcat_value(enz_id, kcat_dict)