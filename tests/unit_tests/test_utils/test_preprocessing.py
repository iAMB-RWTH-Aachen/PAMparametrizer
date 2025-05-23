import pandas as pd
import numpy as np
import pytest

from Modules.utils.preprocessing import (assign_missing_gprs,
                                         map_kcat_values_to_reaction_protein_association)

@pytest.mark.parametrize(
    "input_row,use_ec,expected_gene,expected_gpr,expected_enzyme_id",
    [
        # Case 1: Missing GPR, use_ec = False
        ({'GPR': '', 'rxn_id': 'R1', 'EC': np.nan, 'gene': np.nan, 'enzyme_id': np.nan},
         False, 'Gene_R1', 'Gene_R1', 'Enzyme_R1'),

        # Case 2: Missing GPR, use_ec = True, with EC number
        ({'GPR': '', 'rxn_id': 'R2', 'EC': '1.1.1.1', 'gene': np.nan, 'enzyme_id': np.nan},
         True, 'Gene_R2_1.1.1.1', 'Gene_R2_1.1.1.1', 'Enzyme_R2_1.1.1.1'),

        # Case 3: Missing GPR, use_ec = True, but EC is NaN → should NOT assign
        ({'GPR': '', 'rxn_id': 'R3', 'EC': np.nan, 'gene': np.nan, 'enzyme_id': np.nan},
         True, np.nan, '', np.nan),

        # Case 4: GPR already assigned → should not overwrite
        ({'GPR': 'existing_GPR', 'rxn_id': 'R4', 'EC': '2.2.2.2', 'gene': 'existing_gene', 'enzyme_id': 'existing_enzyme'},
         True, 'existing_gene', 'existing_GPR', 'existing_enzyme'),

        # Case 5: GPR is NaN but use_ec = False → fallback to rxn_id
        ({'GPR': np.nan, 'rxn_id': 'R5', 'EC': '3.3.3.3', 'gene': np.nan, 'enzyme_id': np.nan},
         False, 'Gene_R5', 'Gene_R5', 'Enzyme_R5'),
    ]
)
def test_assign_missing_gprs(input_row, use_ec, expected_gene, expected_gpr, expected_enzyme_id):
    df = pd.DataFrame([input_row])
    result_df = assign_missing_gprs(df, use_ec=use_ec)

    # Check expected values
    assert result_df.loc[0, 'gene'] == expected_gene or (pd.isna(expected_gene) and pd.isna(result_df.loc[0, 'gene']))
    assert result_df.loc[0, 'GPR'] == expected_gpr or (pd.isna(expected_gpr) and pd.isna(result_df.loc[0, 'GPR']))
    assert result_df.loc[0, 'enzyme_id'] == expected_enzyme_id or (
        pd.isna(expected_enzyme_id) and pd.isna(result_df.loc[0, 'enzyme_id']))

@pytest.mark.parametrize("id_mapper, gotenzymes_df, expected", [
    # Case 1: Match by gene ID
    (
        pd.DataFrame({
            'rxn_id': ['R1'],
            'locus_tag': ['gene1'],
            'kegg.reaction': ['RXN1'],
            'EC': ['1.1.1.1'],
            'reversible': [True],
            'Reactants': [['A', 'B']],
            'Products': [['D']]
        }),
        pd.DataFrame({
            'gene': ['gene1'],
            'reaction_id': ['RXN1'],
            'ec_number': ['1.1.1.1'],
            'kcat_values': [10.1]
        }),
        pd.DataFrame({
            'rxn_id': ['R1'],
            'gene': ['gene1'],
            'kcat_values': [10.1]
        })
    ),
    # Case 2: No gene match, EC match only
    (
        pd.DataFrame({
            'rxn_id': ['R2'],
            'locus_tag': ['geneX'],  # Not in GotEnzymes
            'kegg.reaction': ['RXN2'],
            'EC': ['2.2.2.2'],
            'reversible': [True],
            'Reactants': [['A', 'B']],
            'Products': [['D']]
        }),
        pd.DataFrame({
            'gene': ['other_gene'],
            'reaction_id': ['RXN2'],
            'ec_number': ['2.2.2.2'],
            'kcat_values': [22.2]
        }),
        pd.DataFrame({
            'rxn_id': ['R2'],
            'gene': ['geneX'],
            'kcat_values': [22.2]
        })
    ),
    # Case 3: No match at all
    (
        pd.DataFrame({
            'rxn_id': ['R3'],
            'locus_tag': ['geneZ'],
            'kegg.reaction': ['RXN3'],
            'EC': ['3.3.3.3'],
            'reversible': [True],
            'Reactants': [['A', 'B']],
            'Products': [['D']]
        }),
        pd.DataFrame({
            'gene': ['unrelated'],
            'reaction_id': ['OTHER'],
            'ec_number': ['9.9.9.9'],
            'kcat_values': [99.9]
        }),
        pd.DataFrame({
            'rxn_id': ['R3'],
            'gene': ['geneZ'],
            'kcat_values': [np.nan]
        })
    )
])
def test_map_kcat_values_to_reaction_protein_association(id_mapper, gotenzymes_df, expected):
    result = map_kcat_values_to_reaction_protein_association(id_mapper, gotenzymes_df)
    result = result[['rxn_id', 'gene', 'kcat_values']].reset_index(drop=True)
    expected = expected.reset_index(drop=True)
    pd.testing.assert_frame_equal(result, expected)

def test_map_kcat_values_to_reaction_protein_association_preserves_model_columns():
    # Arrange
    id_mapper = pd.DataFrame({
        'rxn_id': ['R1', 'R2'],
        'locus_tag': ['gene1', 'geneX'],
        'kegg.reaction': ['RXN1', 'RXN2'],
        'EC': ['1.1.1.1', '2.2.2.2'],
        'reversible': [True, False],
        'Reactants': [['A', 'B'], ['C']],
        'Products': [['D'], ['E']]
    })

    # GotEnzymes: one row matches by gene, the other only by EC
    gotenzymes_df = pd.DataFrame({
        'gene': ['gene1', 'not_used'],
        'reaction_id': ['RXN1', 'RXN2'],
        'ec_number': ['1.1.1.1', '2.2.2.2'],
        'kcat_values': [10.0, 22.2]
    })

    # Act
    merged = map_kcat_values_to_reaction_protein_association(id_mapper, gotenzymes_df)

    # Assert all required columns are preserved and non-null
    for col in ['rxn_id', 'gene', 'reversible', 'Reactants', 'Products', 'kcat_values']:
        assert col in merged.columns, f"Missing column: {col}"
        assert merged[col].isnull().sum() == 0, f"Null values found in column: {col}"