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

def test_if_parallelized_y_intercept_update_is_saved_in_final_diagnostic_file():
    sut = PAMParametrizerMock()
    sut.result_diagnostics_file = sut.result_diagnostics_file.split('.')[0] + f'test.xlsx'
    ues_config = {'slope': 1e-2, 'intercept': 1e-2}
    sut.validation_data.R1.sector_configs['UnusedEnzymeSector'] = ues_config

    # Act
    sut.optimize_sector_yintercept()
    sut.save_final_diagnostics()

    #Assert
    sector_df = pd.read_excel(sut.result_diagnostics_file, 'sector_parameters')
    assert all([ues_config[k] != v for k,v in sut.validation_data.R1.sector_configs['UnusedEnzymeSector'].items()])
    assert os.path.exists(sut.result_diagnostics_file)
    assert len(sector_df)==2 # Translational Protein and Unused Enzymes sector
    assert sector_df[sector_df.sector_id == 'UnusedEnzymeSector']['slope'].iloc[0] != ues_config['slope']
    assert sector_df[sector_df.sector_id == 'UnusedEnzymeSector']['intercept'].iloc[0] != ues_config['intercept']

    [os.remove(file) for file in [sut.result_diagnostics_file]]

def test_if_y_intercept_update_is_saved_in_final_diagnostic_file():
    sut = PAMParametrizerMock()
    sut.result_diagnostics_file = sut.result_diagnostics_file.split('.')[0] + f'test.xlsx'
    ues_config = {'slope': 1e-2, 'intercept': 1e-2}
    sut.validation_data.R1.sector_configs['UnusedEnzymeSector'] = ues_config

    # Act
    sut. _optimize_sector_yintercept_for_validation_data(vd=sut.validation_data.R1,
                                                         sector_id='UnusedEnzymeSector',
                                                         throw_warning=False)
    sut.save_final_diagnostics()

    #Assert
    sector_df = pd.read_excel(sut.result_diagnostics_file, 'sector_parameters')
    assert all([ues_config[k] != v for k,v in sut.validation_data.R1.sector_configs['UnusedEnzymeSector'].items()])
    assert os.path.exists(sut.result_diagnostics_file)
    assert len(sector_df)==2 # Translational Protein and Unused Enzymes sector
    assert sector_df[sector_df.sector_id == 'UnusedEnzymeSector']['slope'].iloc[0] != ues_config['slope']
    assert sector_df[sector_df.sector_id == 'UnusedEnzymeSector']['intercept'].iloc[0] != ues_config['intercept']

    [os.remove(file) for file in [sut.result_diagnostics_file]]