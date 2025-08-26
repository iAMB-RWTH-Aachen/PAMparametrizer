import pytest
import pandas as pd
from typing import Callable

from Modules.PAM_parametrizer.PAM_data_classes import ParametrizationResults, SectorConfig, ValidationData, FluxResults, KcatConstraintConfigTable
from dataclasses import MISSING, fields

def is_primitive(val):
    return isinstance(val, (int, float, str, bool, type(None)))

@pytest.fixture
def kcatconfig_df():
    return pd.DataFrame({
        "enzyme_id": ["E1", "E2"],
        "reaction_id": ["R1", "R2"],
        "direction": ["f", "b"],
        "min_kcat": [0, 10],
        "max_kcat": [100, 200],
    })


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


def test_kcatconstraintconfig_initializes_from_df(kcatconfig_df):
    cfg = KcatConstraintConfigTable(kcatconfig_df)
    assert isinstance(cfg.df, pd.DataFrame)
    assert ("E1", "R1", "f") in cfg.df.index
    assert cfg.df.loc[("E2", "R2", "b"), "max_kcat"] == 200


def test_kcatconstraintconfig_raises_error_when_missing_columns(kcatconfig_df):
    bad_df = kcatconfig_df.drop(columns=["max_kcat"])
    with pytest.raises(ValueError, match="Missing required columns"):
        KcatConstraintConfigTable(bad_df)

@pytest.mark.parametrize('direction, value, errormessage',[
    ("min", -1, "min_kcat must be >= 0"),
    ("max", 0, "max_kcat must be strictly greater")
])
def test_kcatconstraintconfig_raises_error_with_invalid_kcat(kcatconfig_df, direction, value, errormessage):
    kcatconfig_df.loc[0, f"{direction}_kcat"] = value
    with pytest.raises(ValueError, match=errormessage):
        KcatConstraintConfigTable(kcatconfig_df)


def test_kcatconstraintconfig_gets_valid_value(kcatconfig_df):
    cfg = KcatConstraintConfigTable(kcatconfig_df)
    result = cfg.get("E1", "R1", "f")
    assert result == {"min_kcat": 0, "max_kcat": 100}


def test_kcatconstraintconfig_does_not_get_invalid_value(kcatconfig_df):
    cfg = KcatConstraintConfigTable(kcatconfig_df)
    with pytest.raises(KeyError, match="No kcat constraint"):
        cfg.get("E3", "R3", "f")


def test_kcatconstraintconfig_add_new_constraint(kcatconfig_df):
    cfg = KcatConstraintConfigTable(kcatconfig_df)
    cfg.add("E3", "R3", "f", min_kcat=5, max_kcat=500)
    result = cfg.get("E3", "R3", "f")
    assert result == {"min_kcat": 5, "max_kcat": 500}


def test_kcatconstraintconfig_has_constraint(kcatconfig_df):
    cfg = KcatConstraintConfigTable(kcatconfig_df)
    assert cfg.has_constraint("E1", "R1", "f")
    assert not cfg.has_constraint("E9", "R9", "f")