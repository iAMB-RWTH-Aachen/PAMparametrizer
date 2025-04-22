import pandas as pd
import os
import pytest
from typing import List

from cobra import Reaction
from PAModelpy.utils.pam_generation import parse_reaction2protein

from Modules.utils.pam_generation import _extract_reaction_id_from_catalytic_reaction_id
from Scripts.i2_parametrization.pam_parametrizer_iJN1463 import set_up_pamparametrizer
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
                                           nmbr_kcats_to_pick=len(gpr_relations_to_add))
    enzymes_to_evaluate = sut._parse_enzymes_to_evaluate(esc_topn_df)

    # Assert
    rxn_in_e2evaluate = []
    for rxn in rxns2test:
        #each enzyme can be associated to multiple reactions, and only one association is required per reaction
        rxn_in_e2evaluate += [any([any([rxn in rdict['reaction'] for rdict in edict])
                                   for edict in enzymes_to_evaluate.values()])]
    assert all(rxn_in_e2evaluate)

@pytest.mark.parametrize('sheet_name, rxns2test', [
    ('pseudomonasids', ['ACLS','PEAMNO2pp']),
    ('cglutanicumids', ['4HOXPACt2pp', '2MAHMP'])
])
def test_enzymes_to_evaluate_adds_right_number_of_gprs(sheet_name:str, rxns2test: List[str]):
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
                                           nmbr_kcats_to_pick=len(gpr_relations_to_add))
    enzymes_to_evaluate = sut._parse_enzymes_to_evaluate(esc_topn_df)
    #the number of gprs per enzyme in the enzymes to evaluate depend on how many made it to the topn best results
    # Filter gpr_relations_to_add to only keep valid (rxn, enzyme) pairs that exist in esc_topn_df
    esc_topn_df['reaction_parsed'] = esc_topn_df.rxn_id.apply(
            _extract_reaction_id_from_catalytic_reaction_id
        )
    filtered_gpr = gpr_relations_to_add[
        gpr_relations_to_add[['rxn_id', 'enzyme_id']].apply(tuple, axis=1).isin(
            esc_topn_df[['reaction_parsed', 'enzyme_id']].apply(tuple, axis=1)
        )
    ].drop_duplicates(['rxn_id', 'enzyme_id'], keep ='first')
    # Count the number of GPR relations per enzyme
    num_gpr_per_enzyme = filtered_gpr.groupby('enzyme_id').size().to_dict()


    # Assert
    assert all([len(gpr_list) == num_gpr_per_enzyme[enz]
                for enz, gpr_list in enzymes_to_evaluate.items()
                ])

# def test_iJN1463_parametrization_when_error_is_converging():
#     # Arrange
#     sut = set_up_pamparametrizer(-10,-9)
#     sut.hyperparameters.number_of_kcats_to_mutate = len(sut._pamodel.enzymes)
#     #use dummy ESC results
#     sut._pamodel.optimize()
#     sut.parametrization_results.add_enzyme_sensitivity_coefficients(
#         sut._pamodel.enzyme_sensitivity_coefficients,
#         bin_id='',
#         substrate_uptake_rate=-10)
#
#     # Act
#     enzymes_to_evaluate = sut._determine_enzymes_to_evaluate_for_all_bins(
#         nmbr_kcats_to_pick=sut.hyperparameters.number_of_kcats_to_mutate)

