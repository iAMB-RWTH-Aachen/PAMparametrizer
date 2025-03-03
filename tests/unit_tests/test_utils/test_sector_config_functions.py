import pytest
import numpy as np
import os

from Modules.utils.sector_config_functions import run_simulations
from Modules.utils.pam_generation import set_up_ecoli_pam

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
        assert all([solf == pytest.approx(testf, rel=1e-3) for solf, testf in zip(sol.fluxes, fluxes)])


def test_if_sector_config_run_simulations_gives_same_results_as_running_simulations_no_sectors():
    # Arrange
    pam = set_up_ecoli_pam(total_protein=False, active_enzymes = False, unused_enzymes = False,
                          sensitivity = False)
    substrate_rates = np.arange(-11,-5,1)
    substr_upt_id = 'EX_glc__D_e'

    # Act
    fluxes_to_test = run_simulations(pam, substrate_rates, substr_upt_id)

    # Assert
    for fluxes in fluxes_to_test:
        pam.change_reaction_bounds(substr_upt_id, fluxes[substr_upt_id], 0)
        sol = pam.optimize()
        assert pam.objective.value > 0
        assert all([solf == pytest.approx(testf, rel=1e-3) for solf, testf in zip(sol.fluxes, fluxes)])