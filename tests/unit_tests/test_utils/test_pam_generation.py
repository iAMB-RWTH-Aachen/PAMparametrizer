import os
import pandas as pd
import tempfile
import pytest

from Scripts.pam_generation import setup_toy_pam
from Modules.utils.pam_generation import (create_new_aes_parameter_file,
                                          _extract_reaction_id_from_catalytic_reaction_id,
                                          create_pamodel_from_diagnostics_file)


def test_create_pamodel_from_diagnostics_file_changes_sector():
    #Arrange
    diagnostics_file = os.path.join('tests', 'data', 'diagnostics_file_for_toy.xlsx')
    sheet_name = 'case_01'
    model = setup_toy_pam()
    substrate_uptake_id = 'R1'

    #Act
    create_pamodel_from_diagnostics_file(file_path = diagnostics_file,
                                         model = model,
                                         sheet_name = sheet_name,
                                         substrate_uptake_id = substrate_uptake_id
                                         )

    #Assert
    assert model.sectors.TranslationalProteinSector.slope == 0.000001*1e3 #1e3 to convert g to mg (units of protein abundance)
    assert model.sectors.TranslationalProteinSector.intercept == 0.000001*1e3


def test_if_create_aes_parameter_file_generates_file():
    # Arrange/Act
    temp_dir, files = create_new_aes_file_arrange_act()

    # Assert
    assert os.path.isfile(files['output_file'])
    # cleanup temporary files
    temp_dir.cleanup()

def test_if_create_aes_parameter_file_corectly_sets_up_file():
    # Arrange/Act
    temp_dir, files = create_new_aes_file_arrange_act()
    diagnostics_file = pd.read_excel(files['diagnostics_file'],
                                     sheet_name='case_01')

    # Act
    updated_aes_data = pd.read_excel(files['output_file'],
                                 sheet_name='ActiveEnzymes')

    # Check if kcat values have been updated
    for index,row in diagnostics_file.drop_duplicates(['enzyme_id', 'rxn_id', 'direction'],
                                                      keep = 'last').iterrows():
        rxn_id = _extract_reaction_id_from_catalytic_reaction_id(row['rxn_id'])

        match = updated_aes_data[
            (updated_aes_data["rxn_id"] == rxn_id) &
            (updated_aes_data["enzyme_id"] == row["enzyme_id"]) &
            (updated_aes_data["direction"] == row["direction"])
            ]

        assert not match.empty, f"Missing update for {rxn_id} {row['enzyme_id']}"
        assert match.iloc[0]["kcat_values"] == row[
            "kcat[s-1]"], f"Incorrect kcat update for {rxn_id} {row['enzyme_id']}"

    #cleanup temporary files
    temp_dir.cleanup()

@pytest.mark.parametrize('input_rxn, expected_output', [
    ('PEAMNO2pp', 'PEAMNO2pp' ),
    ('CE_PEAMNO2pp_Enzyme_PP_3462_P0A181_Q88H99_Q88HA1', 'PEAMNO2pp'),
    ('CE_PEAMNO2pp', 'PEAMNO2pp')
])
def test_extract_reaction_id_from_catalytic_reaction_id_gives_right_output(input_rxn:str,
                                                                           expected_output:str):
    # Act
    rxn_id = _extract_reaction_id_from_catalytic_reaction_id(input_rxn)

    assert rxn_id == expected_output


########################################################################################################################
# HELPER FUNCTIONS
########################################################################################################################

def create_new_aes_file_arrange_act():
    # Arrange
    old_param_file = os.path.join('tests', 'data', 'proteinAllocationModel_toymodel.xlsx')
    diagnostics_file = os.path.join('tests', 'data', 'diagnostics_file_for_toy.xlsx')
    # Set up temporary output files for testing.
    temp_dir = tempfile.TemporaryDirectory()
    output_file = os.path.join(temp_dir.name, "new_aes.xlsx")
    files = {'old_param_file': old_param_file,
            'diagnostics_file': diagnostics_file,
            'output_file': output_file}

    # Act
    create_new_aes_parameter_file(old_param_file, diagnostics_file, output_file,
                                  diagnostics_sheet_name = 'case_01')

    return temp_dir, files