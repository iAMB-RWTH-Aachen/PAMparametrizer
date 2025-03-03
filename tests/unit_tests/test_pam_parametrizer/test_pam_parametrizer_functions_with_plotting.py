import os
from tests.pam_parametrizer_mock import PAMParametrizerMock


def test_pam_parametrizer_plots_validation_data():
    # Arrange
    sut = PAMParametrizerMock()

    # Act
    fig, axs = sut.plot_valid_data()


def test_pam_parametrizer_plots_progress():
    # Arrange
    sut = PAMParametrizerMock()
    fig, axs = sut.plot_valid_data()
    # Act
    fig = sut.plot_simulation(fig, axs)


def test_pam_parametrizer_runs_full_workflow_with_bins():
    # Arrange
    sut = PAMParametrizerMock()
    sut.min_substrate_uptake_rate = 0.07
    sut.max_substrate_uptake_rate = 0.09
    sut.hyperparameters.threshold_error = 1

    # Act
    sut.run()
    #remove result files (incl figure)
    # filename_extension = f'final_run_{sut.iteration}'
    # full_file_path = os.path.join('Results', '2_parametrization',
    #                               sut.hyperparameters.genetic_algorithm_filename_base + filename_extension)
    # [os.remove(full_file_path + file_type) for file_type in ['.json', '.xlsx', '.pickle']]
    os.remove(sut.result_figure_file)

    # Assert
    # if it runs all is fine, only testing functionality
    assert True

def test_pam_parametrizer_runs_full_workflow_with_bins_before_iterations():
    # Arrange
    sut = PAMParametrizerMock()
    sut.min_substrate_uptake_rate = 0.07
    sut.max_substrate_uptake_rate = 0.09
    sut.hyperparameters.threshold_error = 1

    # Act
    sut.run(binned = 'before')

    # Assert
    # if it runs all is fine, only testing functionality
    assert True
def test_pam_parametrizer_runs_full_workflow_without_bins():
    # Arrange
    sut = PAMParametrizerMock()
    sut.min_substrate_uptake_rate = 0.07
    sut.max_substrate_uptake_rate = 0.09
    sut.hyperparameters.threshold_error = 1

    # Act
    sut.run(binned = 'False')
    #remove result files (incl figure)
    os.remove(sut.result_figure_file)

    # Assert
    # if it runs all is fine, only testing functionality
    assert True

