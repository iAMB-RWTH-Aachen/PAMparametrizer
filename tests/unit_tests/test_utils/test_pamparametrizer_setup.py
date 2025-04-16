import pytest
import pandas as pd
import os
from typing import List
from Modules.model_setup.sector_config import set_up_sector_config
from Modules.model_setup.sector_classes import SectorConfig
from Modules.utils.sector_config_functions import DictList

@pytest.mark('sector_ids, sector_sheet_names',[
    (['UnusedEnzymeSector'], ['UnusedEnzyme']),
    (['TranslationalProteinSector'], ['Translational']),
    (['UnusedEnzymeSector', 'TranslationalProteinSector'], ['UnusedEnzyme', 'Translational'])
])
def test_set_up_sector_config_returns_expected_structure(sector_ids: List[str],
                                                         sector_sheet_names: List[str]
                                                         ):
    # Arrange
    pam_info_file_path = os.path.join(
        'tests',
        'data',
        'proteinAllocationModel_toymodel.xlsx'
    )
    sectors_info = {}
    for sheet_name, sector_id in zip(sector_sheet_names, sector_ids):
        sectors_info[sector_id]= pd.read_excel(pam_info_file_path,
                                               sheet_name = sheet_name
                                               ).set_index('Parameter')[['Value_for_growth']]

    # Apply
    config = set_up_sector_config(pam_info_file = pam_info_file_path,
                                  sectors_not_related_to_growth = sector_ids
                                  )

    # Assert
    assert isinstance(config, DictList)
    for sector_id in sector_ids:
        assert sector_id in config
        assert isinstance(config[sector_id], SectorConfig)

        sector_info = sectors_info[sector_id]
        assert sector_info.index.str.endswith('_mu') == config[sector_id].slope_vs_mu
        assert sector_info.index.str.endswith('_0') == config[sector_id].intercept_vs_mu
        assert sector_info.substrate_range == config[sector_id].substrate_range
