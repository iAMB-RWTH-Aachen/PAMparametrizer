import pytest
import pandas as pd
from typing import Callable

from Modules.PAM_parametrizer.PAM_data_classes import ParametrizationResults, SectorConfig, ValidationData, FluxResults
import builtins
from dataclasses import MISSING, fields

def is_primitive(val):
    return isinstance(val, (int, float, str, bool, type(None)))

@pytest.mark.parametrize('object, kwargs',[
    (ParametrizationResults, {'substrate_uptake_reactions':['EX_glc__D_e']}),
    # (SectorConfig, {'slope':1, 'intercept':1, 'sectorname': 'TranslationalProteinSector', 'substrate_range':[0,1]}),
    (ValidationData, {'valid_data':pd.DataFrame(),'id': 'EX_glc__D_e',  'validation_range': []}),
    (FluxResults, {'id': 'EX_glc__D_e'})
])
def test_no_shared_mutable_fields_between_data_object_instances(object:Callable, kwargs):
    obj1 = object(**kwargs)
    obj2 = object(**kwargs)
    shared_fields = []
    for f in fields(obj1):
        val1 = getattr(obj1, f.name)
        val2 = getattr(obj2, f.name)

        # --- SKIP RULES ---
        if str(f.name) in kwargs: continue

        if val1 is None and val2 is None:
            continue

        if is_primitive(val1) and is_primitive(val2):
            continue

        # If field has no default_factory and default is a primitive
        if f.default is not MISSING and is_primitive(f.default):
            continue

        # If still same id → suspicious
        if id(val1) == id(val2):
            shared_fields.append(f.name)

    assert not shared_fields, f"Shared fields detected: {shared_fields}"