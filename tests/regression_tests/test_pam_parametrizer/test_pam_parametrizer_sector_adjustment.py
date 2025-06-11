import os
import pandas as pd
import pytest
from time import time

from tests.pam_parametrizer_mock import PAMParametrizerMockEcoli

def test_if_unused_protein_sector_is_configured_correctly_by_pamparametrizer():
    # Arrange
    reference_unused_enzyme_sector = pd.read_excel(
        os.path.join('tests', 'data', 'proteinAllocationModel_iML1515_EnzymaticData_dummy.xlsx'),
        sheet_name = 'UnusedEnzyme'
    ).set_index('Parameter')

    # Apply
    sut = PAMParametrizerMockEcoli()
    unused_enzymes_sut = sut.validation_data.EX_glc__D_e.sector_configs['UnusedEnzymeSector']

    # Assert
    assert reference_unused_enzyme_sector.at['ups_0', 'Value'] == pytest.approx(unused_enzymes_sut['intercept'],
                                                                                abs = 1e-3)
    assert reference_unused_enzyme_sector.at['ups_mu', 'Value'] == pytest.approx(unused_enzymes_sut['slope'],
                                                                                 abs = 1e-3)

def test_if_yintercept_UE_is_correctly_reconfigured_by_pamparametrizer():
    #Arrange
    sut = PAMParametrizerMockEcoli()

    # Act
    sut.optimize_sector_yintercept(sector_id='UnusedEnzymeSector')

    # Assert
    for vd in sut.validation_data:
        params = vd.sector_configs['UnusedEnzymeSector']
        assert params['intercept'] >= sut.minimal_unused_enzymes * sut._pamodel.total_protein_fraction
        assert isinstance(params['slope'], float)

def test_if_yintercept_UE_is_fast_enough():
    #Arrange
    sut = PAMParametrizerMockEcoli()
    start = time()
    # Act
    sut.optimize_sector_yintercept(sector_id='UnusedEnzymeSector')
    assert time()-start < 15*60 #less than 15 min

