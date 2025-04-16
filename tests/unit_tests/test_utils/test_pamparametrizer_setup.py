import ast
import pytest
import pandas as pd
import os
from typing import List

from Modules.utils.pamparametrizer_setup import set_up_sector_config
from Modules.PAM_parametrizer import SectorConfig

@pytest.mark.parametrize('sector_ids, sector_sheet_names',[
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
    assert isinstance(config, dict)
    for sector_id in sector_ids:
        assert sector_id in config
        sector_info = sectors_info[sector_id]
        substrate_range_expected = ast.literal_eval(sector_info.loc['substrate_range', 'Value_for_growth'])
        assert sector_info.filter(like='_mu', axis=0)['Value_for_growth'].iloc[0] == config[sector_id]['slope']
        assert sector_info.filter(like='_0', axis=0)['Value_for_growth'].iloc[0] == config[sector_id]['intercept']
        assert substrate_range_expected == config[sector_id]['substrate_range']
