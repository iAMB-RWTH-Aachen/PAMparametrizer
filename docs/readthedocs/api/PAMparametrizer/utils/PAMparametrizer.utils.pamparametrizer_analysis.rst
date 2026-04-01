:py:mod:`PAMparametrizer.utils.pamparametrizer_analysis`
========================================================

.. py:module:: PAMparametrizer.utils.pamparametrizer_analysis

.. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`set_up_pam_parametrizer_and_get_substrate_uptake_rates <PAMparametrizer.utils.pamparametrizer_analysis.set_up_pam_parametrizer_and_get_substrate_uptake_rates>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.set_up_pam_parametrizer_and_get_substrate_uptake_rates
          :summary:
   * - :py:obj:`get_results_from_simulations <PAMparametrizer.utils.pamparametrizer_analysis.get_results_from_simulations>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.get_results_from_simulations
          :summary:
   * - :py:obj:`get_results_from_simulations_fixed_mu <PAMparametrizer.utils.pamparametrizer_analysis.get_results_from_simulations_fixed_mu>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.get_results_from_simulations_fixed_mu
          :summary:
   * - :py:obj:`_set_up_pamodel_for_simulations <PAMparametrizer.utils.pamparametrizer_analysis._set_up_pamodel_for_simulations>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis._set_up_pamodel_for_simulations
          :summary:
   * - :py:obj:`_set_up_solution_info <PAMparametrizer.utils.pamparametrizer_analysis._set_up_solution_info>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis._set_up_solution_info
          :summary:
   * - :py:obj:`save_fluxes <PAMparametrizer.utils.pamparametrizer_analysis.save_fluxes>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.save_fluxes
          :summary:
   * - :py:obj:`save_proteins <PAMparametrizer.utils.pamparametrizer_analysis.save_proteins>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.save_proteins
          :summary:
   * - :py:obj:`calculate_error_for_reactions <PAMparametrizer.utils.pamparametrizer_analysis.calculate_error_for_reactions>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.calculate_error_for_reactions
          :summary:
   * - :py:obj:`calculate_r_squared_for_reaction <PAMparametrizer.utils.pamparametrizer_analysis.calculate_r_squared_for_reaction>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.calculate_r_squared_for_reaction
          :summary:
   * - :py:obj:`calculate_difference_simulation_experiment <PAMparametrizer.utils.pamparametrizer_analysis.calculate_difference_simulation_experiment>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.calculate_difference_simulation_experiment
          :summary:
   * - :py:obj:`calculate_kcat_differences <PAMparametrizer.utils.pamparametrizer_analysis.calculate_kcat_differences>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.calculate_kcat_differences
          :summary:
   * - :py:obj:`calculate_relative_change <PAMparametrizer.utils.pamparametrizer_analysis.calculate_relative_change>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.calculate_relative_change
          :summary:
   * - :py:obj:`get_previous_kcat_values <PAMparametrizer.utils.pamparametrizer_analysis.get_previous_kcat_values>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.get_previous_kcat_values
          :summary:
   * - :py:obj:`convert_peptide_to_enzyme_concentrations <PAMparametrizer.utils.pamparametrizer_analysis.convert_peptide_to_enzyme_concentrations>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.convert_peptide_to_enzyme_concentrations
          :summary:
   * - :py:obj:`parse_enzyme_complex_id <PAMparametrizer.utils.pamparametrizer_analysis.parse_enzyme_complex_id>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.parse_enzyme_complex_id
          :summary:
   * - :py:obj:`normalize_simulated_protein_concentrations <PAMparametrizer.utils.pamparametrizer_analysis.normalize_simulated_protein_concentrations>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.normalize_simulated_protein_concentrations
          :summary:
   * - :py:obj:`get_clusters_from_clustermap <PAMparametrizer.utils.pamparametrizer_analysis.get_clusters_from_clustermap>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.get_clusters_from_clustermap
          :summary:
   * - :py:obj:`select_clustered_rows_by_variation <PAMparametrizer.utils.pamparametrizer_analysis.select_clustered_rows_by_variation>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.select_clustered_rows_by_variation
          :summary:
   * - :py:obj:`row_wise_zscore_normalization <PAMparametrizer.utils.pamparametrizer_analysis.row_wise_zscore_normalization>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.row_wise_zscore_normalization
          :summary:
   * - :py:obj:`_compute_variation_metrics_per_row <PAMparametrizer.utils.pamparametrizer_analysis._compute_variation_metrics_per_row>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis._compute_variation_metrics_per_row
          :summary:
   * - :py:obj:`plot_histogram_logspace <PAMparametrizer.utils.pamparametrizer_analysis.plot_histogram_logspace>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.plot_histogram_logspace
          :summary:
   * - :py:obj:`plot_PCA_graph <PAMparametrizer.utils.pamparametrizer_analysis.plot_PCA_graph>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.plot_PCA_graph
          :summary:

API
~~~

.. py:function:: set_up_pam_parametrizer_and_get_substrate_uptake_rates(set_up_parametrizer: typing.Callable, parametrizer_kwargs: typing.Dict = {'max_substrate_uptake_rate': -0.1}, substrate_uptake_id: str = 'EX_glc__D_e') -> typing.Tuple
   :canonical: PAMparametrizer.utils.pamparametrizer_analysis.set_up_pam_parametrizer_and_get_substrate_uptake_rates

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.set_up_pam_parametrizer_and_get_substrate_uptake_rates

.. py:function:: get_results_from_simulations(model: typing.Union[PAModelpy.PAModel, Model], substrate_rates: typing.Union[typing.Iterable[float], typing.Iterable[typing.Iterable[float]]], substrate_ids: typing.Union[str, list[str]] = ['EX_glc__D_e'], fluxes_to_save: list[str] = None, proteins_to_save: list[str] = None, sectors_config: typing.Union[dict[str, PAMparametrizer.utils.sector_config_functions.SectorParameterDict], bool] = True) -> dict[str, pandas.DataFrame]
   :canonical: PAMparametrizer.utils.pamparametrizer_analysis.get_results_from_simulations

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.get_results_from_simulations

.. py:function:: get_results_from_simulations_fixed_mu(pamodel: PAModelpy.PAModel, growth_rates: typing.Iterable[float], substrate_id: str = 'EX_glc__D_e', fluxes_to_save: list[str] = None, proteins_to_save: list[str] = None, sector_config=True, method_ids: list[str] = None) -> dict[str, pandas.DataFrame]
   :canonical: PAMparametrizer.utils.pamparametrizer_analysis.get_results_from_simulations_fixed_mu

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.get_results_from_simulations_fixed_mu

.. py:function:: _set_up_pamodel_for_simulations(pamodel: PAModelpy.PAModel, substrate_id: str, sectors_config: typing.Union[bool, typing.Dict[str, PAMparametrizer.utils.sector_config_functions.SectorParameterDict]]) -> None
   :canonical: PAMparametrizer.utils.pamparametrizer_analysis._set_up_pamodel_for_simulations

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis._set_up_pamodel_for_simulations

.. py:function:: _set_up_solution_info(fluxes_to_save: list[str], proteins_to_save: list[str], method_ids: list[str] = None) -> dict[str, pandas.DataFrame]
   :canonical: PAMparametrizer.utils.pamparametrizer_analysis._set_up_solution_info

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis._set_up_solution_info

.. py:function:: save_fluxes(solution: cobra.Solution, pamodel: PAModelpy.PAModel, fluxes_to_save: list[str], substrate_rate: float, flux_df: pandas.DataFrame, substrate_id: str, method: str = None) -> pandas.DataFrame
   :canonical: PAMparametrizer.utils.pamparametrizer_analysis.save_fluxes

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.save_fluxes

.. py:function:: save_proteins(pamodel: PAModelpy.PAModel, proteins_to_save: list[str], substrate_rate: float, protein_df: pandas.DataFrame, substrate_id: str, method: str = None) -> pandas.DataFrame
   :canonical: PAMparametrizer.utils.pamparametrizer_analysis.save_proteins

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.save_proteins

.. py:function:: calculate_error_for_reactions(validation_df: pandas.DataFrame, flux_df: pandas.DataFrame, rxns_to_validate: list) -> float
   :canonical: PAMparametrizer.utils.pamparametrizer_analysis.calculate_error_for_reactions

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.calculate_error_for_reactions

.. py:function:: calculate_r_squared_for_reaction(reaction_id: str, validation_data: pandas.DataFrame, substrate_uptake_id: str, fluxes: pandas.DataFrame) -> float
   :canonical: PAMparametrizer.utils.pamparametrizer_analysis.calculate_r_squared_for_reaction

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.calculate_r_squared_for_reaction

.. py:function:: calculate_difference_simulation_experiment(validation_df: pandas.DataFrame, flux_df: pandas.DataFrame, rxns_to_validate: list[str], substr_rxn: str)
   :canonical: PAMparametrizer.utils.pamparametrizer_analysis.calculate_difference_simulation_experiment

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.calculate_difference_simulation_experiment

.. py:function:: calculate_kcat_differences(df_grouped: pandas.DataFrame, rxn2kcat: dict[str, dict]) -> pandas.DataFrame
   :canonical: PAMparametrizer.utils.pamparametrizer_analysis.calculate_kcat_differences

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.calculate_kcat_differences

.. py:function:: calculate_relative_change(group: pandas.DataFrame, rxn2kcat: dict[str, dict]) -> pandas.DataFrame
   :canonical: PAMparametrizer.utils.pamparametrizer_analysis.calculate_relative_change

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.calculate_relative_change

.. py:function:: get_previous_kcat_values(group: pandas.DataFrame, rxn2kcat: dict[str, dict]) -> float
   :canonical: PAMparametrizer.utils.pamparametrizer_analysis.get_previous_kcat_values

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.get_previous_kcat_values

.. py:function:: convert_peptide_to_enzyme_concentrations(peptide_df: pandas.DataFrame, enzyme_db: pandas.DataFrame, concentration_columns: typing.List[str]) -> pandas.DataFrame
   :canonical: PAMparametrizer.utils.pamparametrizer_analysis.convert_peptide_to_enzyme_concentrations

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.convert_peptide_to_enzyme_concentrations

.. py:function:: parse_enzyme_complex_id(enzyme_id: str)
   :canonical: PAMparametrizer.utils.pamparametrizer_analysis.parse_enzyme_complex_id

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.parse_enzyme_complex_id

.. py:function:: normalize_simulated_protein_concentrations(df: pandas.DataFrame, enzyme_db: pandas.DataFrame, ue_sector: PAModelpy.UnusedEnzymeSector)
   :canonical: PAMparametrizer.utils.pamparametrizer_analysis.normalize_simulated_protein_concentrations

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.normalize_simulated_protein_concentrations

.. py:function:: get_clusters_from_clustermap(clustermap: seaborn.matrix.ClusterGrid, df: pandas.DataFrame, nrow_clusters: int = 10, ncol_clusters: int = 2)
   :canonical: PAMparametrizer.utils.pamparametrizer_analysis.get_clusters_from_clustermap

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.get_clusters_from_clustermap

.. py:function:: select_clustered_rows_by_variation(clustered_data: pandas.DataFrame, column_clusters: dict[str, list], per_cluster: bool = True, select_highest: bool = True, num_rows: int = 10, metric: typing.Literal[MAD, ENT, STD, CV] = 'MAD') -> pandas.DataFrame
   :canonical: PAMparametrizer.utils.pamparametrizer_analysis.select_clustered_rows_by_variation

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.select_clustered_rows_by_variation

.. py:function:: row_wise_zscore_normalization(df: pandas.DataFrame)
   :canonical: PAMparametrizer.utils.pamparametrizer_analysis.row_wise_zscore_normalization

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.row_wise_zscore_normalization

.. py:function:: _compute_variation_metrics_per_row(normalized_data: pandas.DataFrame, column_clusters: dict[str, list]) -> pandas.DataFrame
   :canonical: PAMparametrizer.utils.pamparametrizer_analysis._compute_variation_metrics_per_row

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis._compute_variation_metrics_per_row

.. py:function:: plot_histogram_logspace(ax: matplotlib.pyplot.Axes, data: pandas.DataFrame, color: typing.Union[str, float], label: str, x_label: str, n_bins: int = 50, relative: bool = False, ymax: float = None) -> None
   :canonical: PAMparametrizer.utils.pamparametrizer_analysis.plot_histogram_logspace

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.plot_histogram_logspace

.. py:function:: plot_PCA_graph(input_df, columns_to_analyse: list[str], hue: str, values: str)
   :canonical: PAMparametrizer.utils.pamparametrizer_analysis.plot_PCA_graph

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_analysis.plot_PCA_graph
