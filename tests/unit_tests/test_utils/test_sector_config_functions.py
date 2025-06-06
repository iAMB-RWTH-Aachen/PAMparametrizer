import pytest
import numpy as np
import os

from PAModelpy.utils.pam_generation import set_up_pam
from Modules.utils.sector_config_functions import run_simulations, get_translational_sector_config
from Modules.utils.pam_generation import _extract_reaction_id_from_catalytic_reaction_id
from Scripts.pam_generation import setup_toy_pam

def test_if_sector_config_run_simulations_gives_same_results_as_running_simulations():
    # Arrange
    pam = set_up_ecoli_pam(sensitivity=False)
    substrate_rates = np.arange(-11,-5,1)
    substr_upt_id = 'EX_glc__D_e'

    # Act
    fluxes_to_test = run_simulations(pam, substrate_rates, substr_upt_id)

    # Assert
    for fluxes in fluxes_to_test:
        pam.change_reaction_bounds(substr_upt_id, fluxes[substr_upt_id], 0)
        sol = pam.optimize()

        assert pam.objective.value > 0
        for flux_id, testf in fluxes.items():
            rxn_id = _extract_reaction_id_from_catalytic_reaction_id(flux_id) #to avoid clashed with very similar isozymes
            assert sol[rxn_id] != pytest.approx(fluxes[rxn_id], rel=1e-3)


def test_get_translational_sector_config_returns_sector_correctly():
    #Arrange
    pam = setup_toy_pam(sensitivity=False)
    substrate_id = 'R1'
    substrate_range = np.arange(1e-3, 1e-2, 1e-3)
    expected_translational_sector = {'slope': 0.01*1e-3, 'intercept': 0.01*1e-3}

    #Act
    trans_sector_config = get_translational_sector_config(pam,
                                                          substrate_id,
                                                          substrate_range,
                                                          'R1')
    #Assert
    assert all([trans_sector_config[key] == pytest.approx(value, rel=1e-3) for key, value in expected_translational_sector.items()])


########################################################################################################################
# HELPER FUNCTIONS
########################################################################################################################
def set_up_ecoli_pam(total_protein = 0.258,
                     active_enzymes = True,
                     unused_enzymes = True,
                     sensitivity = False):
    return set_up_pam(
        pam_info_file = os.path.join(
            'tests', 'data', 'proteinAllocationModel_iML1515_EnzymaticData_dummy.xlsx'
        ),
        total_protein = total_protein,
        active_enzymes = True,
        unused_enzymes = True,
        sensitivity = sensitivity)
