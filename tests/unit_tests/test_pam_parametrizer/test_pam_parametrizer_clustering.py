from tests.unit_tests.test_pam_parametrizer.pam_parametrizer_mock import PAMParametrizerMock

# def test_pam_parametrizer_clusters_results_correctly():
#     # Arrange
#     sut = PAMParametrizerMock()
#     parameter_sets = [[1, 2, 3], [1, 4, 3], [7, 8, 9], [2, 4, 6], [8, 7, 5]]
#     fitness_values = [0.9, 0.8, 0.7, 0.6, 0.85]
#     n_clusters = 2
#     expected_cluster_labels = [1, 1, 0,1,0]
#     expected_cluster_fitness = [0.7749999999999999, 0.7666666666666667]
#
#     # Act
#     clusters = sut.cluster_parametrization_results(parameter_sets, fitness_values, n_clusters)
#
#     # Assert
#     print(clusters)
#     # print("Estimated fitness values for each cluster:", clustered_fitness)
    # assert expected_cluster_labels == cluster_labels
    # assert expected_cluster_fitness == clustered_fitness