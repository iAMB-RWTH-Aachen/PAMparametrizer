import pytest
import pandas as pd
from typing import Callable

from Modules.PAM_parametrizer.PAM_data_classes import ParametrizationResults, SectorConfig, ValidationData, FluxResults

@pytest.mark.parametrize('object, kwargs',[
    (ParametrizationResults, {'substrate_uptake_reactions':['EX_glc__D_e']}),
    (SectorConfig, {'slope':1, 'intercept':1, 'sectorname': 'TranslationalProteinSector', 'substrate_range':[0,1]}),
    (ValidationData, {'valid_data':pd.DataFrame(),'id': 'EX_glc__D_e',  'validation_range': []}),
    (FluxResults, {'id': 'EX_glc__D_e'})
])
def test_if_ParametrizationResults_is_initialized_properly(object:Callable, kwargs):
    obj1 = object(**kwargs)
    obj2 = object(**kwargs)
    for attr in obj1.__dataclass_fields__:
        assert id(getattr(obj1, attr)) == id(getattr(obj2, attr)),f"Warning: Field '{attr}' is shared between instances!"