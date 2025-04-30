import os
import pandas as pd
import pytest

from tests.pam_parametrizer_mock import PAMParametrizerMock



# def test_if_diagnostics_are_saved_properly_for_single_run():
#     # Arrange
#     sut = PAMParametrizerMock()
#
#     # Act
#     sut.run()
#     best_indiv_sheet = pd.read_excel(sut.result_diagnostics_file, sheet_name = 'Best_Individuals')
#
#     #Assert
#     assert os.path.exists(sut.result_diagnostics_file)
#     assert len(best_indiv_sheet['run_id'].unique())==3
#     for run, enzyme_df in best_indiv_sheet.groupby('run_id'):
#         assert len(enzyme_df)>=5
#
#     [os.remove(file) for file in [sut.result_diagnostics_file, sut.result_figure_file]]

def test_if_diagnostics_are_saved_properly_for_multiple_runs():
    for iteration in range(2):
        # Arrange
        sut = PAMParametrizerMock()
        sut.result_diagnostics_file = sut.result_diagnostics_file.split('.')[0] + f'test{iteration}.xlsx'

        # Act
        print(f"Iteration {iteration}: ID of parametrization_results = {id(sut.parametrization_results)}")
        sut.run()
        best_indiv_sheet = pd.read_excel(sut.result_diagnostics_file, sheet_name = 'Best_Individuals')

        #Assert
        assert os.path.exists(sut.result_diagnostics_file)
        assert len(best_indiv_sheet['run_id'].unique())==3
        for run, enzyme_df in best_indiv_sheet.groupby('run_id'):
            assert len(enzyme_df)>=5

        [os.remove(file) for file in [sut.result_diagnostics_file, sut.result_figure_file]]