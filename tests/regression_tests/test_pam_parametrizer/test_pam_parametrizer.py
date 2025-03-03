import pandas as pd
import os
import numpy as np
import pytest
from typing import List

from cobra import Reaction
from PAModelpy.utils.pam_generation import parse_reaction2protein
from tests.pam_parametrizer_mock import PAMParametrizerMock

@pytest.mark.parametrize('sheet_name, rxns2test', [
    ('pseudomonasids', ['ACLS','PEAMNO2pp']),
    ('cglutanicumids', ['4HOXPACt2pp', '2MAHMP'])
])
def test_enzymes_to_evaluate_are_parsed_correctly(sheet_name:str, rxns2test: List[str]):
    # Arrange
    sut = PAMParametrizerMock()
    #make sure the information to test is associated with the toy model
    gpr_relations_to_add = pd.read_excel(os.path.join('tests', 'data','mock_enzyme_sensitivities.xlsx'),
                                   sheet_name=sheet_name+'_gpr')
    sut._pamodel.add_reactions([Reaction(rid) for rid in gpr_relations_to_add.rxn_id.drop_duplicates()])
    rxn2protein, protein2gpr = parse_reaction2protein(enzyme_db = gpr_relations_to_add,
                                                      model =sut._pamodel)
    sut._pamodel.add_rxn2protein_to_active_enzymes(rxn2protein, protein2gpr, verbose=True)

    # get dummy enzyme sensitivities
    esc_results_df = pd.read_excel(os.path.join('tests', 'data','mock_enzyme_sensitivities.xlsx'),
                                   sheet_name=sheet_name)

    # Act
    esc_topn_df = sut._select_topn_enzymes(esc_results_df,
                                           nmbr_kcats_to_pick=len(esc_results_df))
    enzymes_to_evaluate = sut._parse_enzymes_to_evaluate(esc_topn_df)

    # Assert
    rxn_in_e2evaluate = []
    for rxn in rxns2test:
        #each enzyme can be associated to multiple reactions, and only one association is required per reaction
        rxn_in_e2evaluate += [any([any([rxn in rdict['reaction'] for rdict in edict])
                                   for edict in enzymes_to_evaluate.values()])]
    assert all(rxn_in_e2evaluate)
