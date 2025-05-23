import pandas as pd
import numpy as np
from Modules.utils.preprocessing import (assign_directionalities_for_kcat_relations,
                                         assign_defaults_for_proteins_without_mapping,
                                         assign_missing_gprs)

DEFAULT_KCAT = 99.9
DEFAULT_MOLMASS = 40000
DEFAULT_PROT_LENGTH = 300


def test_assign_directionalities_for_kcat_relations():
    df = pd.DataFrame({
        'rxn_id': ['R1', 'R1', 'R2', 'R3'],
        'compound': ['C1', 'C2', 'C3', np.nan],
        'Products': [['C1'], ['C4'], ['C3'],['C6']],
        'Reactants': [['C5'], ['C2'], ['C6'],['C1']],
        'reversible': [True, True, False, True],
        'kcat_values': [1.2, 1.5, 2.1, np.nan]
    })

    result = assign_directionalities_for_kcat_relations(df)
    directions_r1 = result[result['rxn_id'] == 'R1']['direction'].values

    # Original 3 rows plus 1 duplicated for the backward direction of R1
    assert len(result) == len(df)+1

    # Ensure direction column is set properly
    assert set(result['direction']).issubset({'f', 'b'})

    # R1 should have both directions
    assert 'b' in directions_r1
    assert 'f' in directions_r1

    # R2 is irreversible, should not be duplicated
    assert result[result['rxn_id'] == 'R2'].shape[0] == 1


def test_assign_defaults_for_proteins_without_mapping(monkeypatch):
    df = pd.DataFrame({
        'kcat_values': [np.nan, 5.0],
        'Mass': [np.nan, 50000],
        'Length': [np.nan, 350],
        'gene': [np.nan, 'geneB'],
        'enzyme_id': [np.nan, 'Enzyme_geneB'],
        'GPR': ['s0001', 'GPR2'],
        'rxn_id': ['R1', 'R2'],
        'EC': ['1.1.1.1', np.nan],
        'compound': ['C1', 'C2'],
        'ec_number': [np.nan, np.nan]
    })

    result = assign_defaults_for_proteins_without_mapping(
        df,
        default_kcat=DEFAULT_KCAT,
        default_molmass=DEFAULT_MOLMASS,
        default_protein_length=DEFAULT_PROT_LENGTH
    )

    # Check default values were assigned
    assert result.loc[0, 'kcat_values'] == DEFAULT_KCAT
    assert result.loc[0, 'molMass'] == DEFAULT_MOLMASS
    assert result.loc[0, 'Length'] == DEFAULT_PROT_LENGTH

    # Check GPR fallback assignment (mocked)
    assert 'gene_from_ec_1.1.1.1' in result.loc[0, 'gene'] or result.loc[0, 'gene'] == 's0001'

    # Check special s0001 handling
    assert result.loc[0, 'enzyme_id'].startswith('Enzyme_s0001R1')

    # Check dropped columns
    assert 'ec_number' not in result.columns
    assert 'compound' not in result.columns

