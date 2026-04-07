:py:mod:`PAMparametrizer.PAM_parametrizer.pam_parametrizer`
===========================================================

.. py:module:: PAMparametrizer.PAM_parametrizer.pam_parametrizer

.. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer
   :allowtitles:

Module Contents
---------------

Classes
~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`PAMParametrizer <PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer>`
     - .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`Reaction2KcatDict <PAMparametrizer.PAM_parametrizer.pam_parametrizer.Reaction2KcatDict>`
     - .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.Reaction2KcatDict
          :summary:
   * - :py:obj:`Enzyme2Reaction2KcatDict <PAMparametrizer.PAM_parametrizer.pam_parametrizer.Enzyme2Reaction2KcatDict>`
     - .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.Enzyme2Reaction2KcatDict
          :summary:

API
~~~

.. py:data:: Reaction2KcatDict
   :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.Reaction2KcatDict
   :value: None

   .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.Reaction2KcatDict

.. py:data:: Enzyme2Reaction2KcatDict
   :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.Enzyme2Reaction2KcatDict
   :value: None

   .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.Enzyme2Reaction2KcatDict

.. py:class:: PAMParametrizer(pamodel: PAModelpy.PAModel.PAModel, validation_data: typing.Union[cobra.DictList[PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData], list, PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData], kcat_configs: typing.Optional[typing.Union[PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable, pandas.DataFrame]] = None, hyperparameters: typing.Optional[PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters] = None, sector_configs: typing.Optional[typing.Union[typing.Dict[str, PAMparametrizer.PAM_parametrizer.PAM_data_classes.SectorConfig], bool]] = None, substrate_uptake_id: str = 'EX_glc__D_e', max_substrate_uptake_rate: typing.Union[float, int] = 0, min_substrate_uptake_rate: typing.Union[float, int] = -11, sensitivity: bool = True, minimal_unused_enzymes: float = 0.37, enzymes_to_evaluate: typing.Optional[list] = None)
   :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer

   .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer

   .. rubric:: Initialization

   .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.__init__

   .. py:attribute:: TRANSLATIONAL_SECTOR_CONFIG
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.TRANSLATIONAL_SECTOR_CONFIG
      :value: 'SectorConfig(...)'

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.TRANSLATIONAL_SECTOR_CONFIG

   .. py:attribute:: TRANSL_SECTOR_INTERCEPT_vs_GLC
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.TRANSL_SECTOR_INTERCEPT_vs_GLC
      :value: 0.046136644909661115

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.TRANSL_SECTOR_INTERCEPT_vs_GLC

   .. py:attribute:: TRANSL_SECTOR_SLOPE_vs_GLC
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.TRANSL_SECTOR_SLOPE_vs_GLC
      :value: None

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.TRANSL_SECTOR_SLOPE_vs_GLC

   .. py:attribute:: MEASURED_PROTEIN_FRACTION
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.MEASURED_PROTEIN_FRACTION
      :value: None

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.MEASURED_PROTEIN_FRACTION

   .. py:attribute:: MAXIMAL_INTERCEPT_UE
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.MAXIMAL_INTERCEPT_UE
      :value: 0.67

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.MAXIMAL_INTERCEPT_UE

   .. py:property:: sector_parameters
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.sector_parameters

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.sector_parameters

   .. py:method:: run(remove_subruns: bool = True, binned: str = 'False') -> None
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.run

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.run

   .. py:method:: calculate_sector_parameters_for_multiple_csources(reset: bool = False)
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.calculate_sector_parameters_for_multiple_csources

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.calculate_sector_parameters_for_multiple_csources

   .. py:method:: perform_iteration_in_bins(start_time: float) -> list
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.perform_iteration_in_bins

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.perform_iteration_in_bins

   .. py:method:: perform_iteration_without_bins(random: bool = False) -> list
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.perform_iteration_without_bins

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.perform_iteration_without_bins

   .. py:method:: optimize_sector_yintercept(sector_id: str = 'UnusedEnzymeSector', throw_warning: bool = False)
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.optimize_sector_yintercept

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.optimize_sector_yintercept

   .. py:method:: _optimize_sector_yintercept_for_validation_data(vd: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData, sector_id: str, throw_warning: bool) -> typing.Tuple[str, PAMparametrizer.PAM_parametrizer.PAM_data_classes.SectorConfig]
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._optimize_sector_yintercept_for_validation_data

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._optimize_sector_yintercept_for_validation_data

   .. py:method:: evaluate_and_save_results_of_iteration(start_time_iteration: float, files_to_remove: list, remove_subruns: bool, fig: matplotlib.pyplot.Figure, axs: matplotlib.pyplot.Axes)
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.evaluate_and_save_results_of_iteration

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.evaluate_and_save_results_of_iteration

   .. py:method:: process_bin(bin_id: str, bin_information: list, substrate_uptake_id: str) -> None
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.process_bin

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.process_bin

   .. py:method:: run_pamodel_simulations_in_bin(substrate_uptake_reaction: str, bin_id: typing.Union[float, int], bin_information: list = None, substrate_uptake_rates: typing.Union[list, numpy.array, pandas.Series] = None) -> None
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.run_pamodel_simulations_in_bin

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.run_pamodel_simulations_in_bin

   .. py:method:: run_pamodel_simulation(substrate_uptake_rate: typing.Union[float, int], bin_id: typing.Union[str, float, int], substrate_uptake_reaction: str) -> None
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.run_pamodel_simulation

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.run_pamodel_simulation

   .. py:method:: save_pamodel_simulation_results(substrate_uptake_rate: typing.Union[float, int], bin_id: typing.Union[str, float, int], substrate_uptake_reaction: str) -> None
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.save_pamodel_simulation_results

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.save_pamodel_simulation_results

   .. py:method:: calculate_error(bin_id: typing.Union[float, int, str], substrate_uptake_reaction: str) -> None
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.calculate_error

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.calculate_error

   .. py:method:: calculate_final_error() -> float
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.calculate_final_error

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.calculate_final_error

   .. py:method:: _calculate_error_for_validation_data(valid_data: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData, bin_id: str = 'final') -> float
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._calculate_error_for_validation_data

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._calculate_error_for_validation_data

   .. py:method:: determine_most_sensitive_enzymes(bin_id: typing.Union[str, float, int], nmbr_kcats_to_pick: int) -> pandas.DataFrame
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.determine_most_sensitive_enzymes

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.determine_most_sensitive_enzymes

   .. py:method:: determine_bin_to_split(esc_topn_df: pandas.DataFrame, bin_id: typing.Union[str, float, int]) -> None
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.determine_bin_to_split

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.determine_bin_to_split

   .. py:method:: run_genetic_algorithm(enzymes_to_evaluate: PAMparametrizer.PAM_parametrizer.pam_parametrizer.Enzyme2Reaction2KcatDict, filename_extension: str) -> None
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.run_genetic_algorithm

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.run_genetic_algorithm

   .. py:method:: restart_genetic_algorithm() -> None
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.restart_genetic_algorithm

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.restart_genetic_algorithm

   .. py:method:: reparametrize_pam(best_individual_kcat_df: pandas.DataFrame = None) -> None
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.reparametrize_pam

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.reparametrize_pam

   .. py:method:: revert_parametrization() -> None
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.revert_parametrization

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.revert_parametrization

   .. py:method:: save_diagnostics(computational_time: float, results_filename: str = None) -> None
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.save_diagnostics

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.save_diagnostics

   .. py:method:: save_final_diagnostics(figure: matplotlib.pyplot.Figure = None) -> None
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.save_final_diagnostics

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.save_final_diagnostics

   .. py:method:: add_new_substrate_source(new_substrate_uptake_id: str, validation_data: pandas.DataFrame, substrate_range: list[typing.Union[int, float]], reactions_to_validate: list)
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.add_new_substrate_source

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.add_new_substrate_source

   .. py:property:: pamodel
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.pamodel
      :type: PAModelpy.PAModel.PAModel

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.pamodel

   .. py:method:: _set_pamodel_no_sensitivities() -> None
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._set_pamodel_no_sensitivities

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._set_pamodel_no_sensitivities

   .. py:method:: _create_result_dirs(result_file_path: str = os.path.join('Results', '2_parametrization')) -> None
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._create_result_dirs

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._create_result_dirs

   .. py:method:: _change_sector_parameters_for_substrate(substrate_uptake_id: str, pamodel: PAModelpy.PAModel.PAModel) -> None
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._change_sector_parameters_for_substrate

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._change_sector_parameters_for_substrate

   .. py:method:: _get_substrate_range_lower_substrate_conc(validation_range: list[typing.Union[int, float]], number_of_steps: int = 5) -> typing.Iterable
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._get_substrate_range_lower_substrate_conc

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._get_substrate_range_lower_substrate_conc

   .. py:method:: _pamodel_is_feasible(substrate_uptake_rate: typing.Optional[float] = None) -> bool
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._pamodel_is_feasible

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._pamodel_is_feasible

   .. py:method:: _bin_substrate_uptake_rates()
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._bin_substrate_uptake_rates

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._bin_substrate_uptake_rates

   .. py:method:: _make_new_bin(bin_id: typing.Union[float, int], bin_range: typing.Union[float, int], substrate_start: typing.Union[float, int])
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._make_new_bin

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._make_new_bin

   .. py:method:: _adjust_binsize(bin_id, bin_range, substrate_start)
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._adjust_binsize

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._adjust_binsize

   .. py:method:: _init_genetic_algorithm(substrate_uptake_rates: dict, enzymes_to_evaluate: PAMparametrizer.PAM_parametrizer.pam_parametrizer.Enzyme2Reaction2KcatDict, sector_configs_per_substrate: typing.Dict[str, PAMparametrizer.utils.sector_config_functions.SectorParameterDict], filename_extension: str) -> PAMparametrizer.genetic_algorithm_parametrization.GAPOGaussian
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._init_genetic_algorithm

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._init_genetic_algorithm

   .. py:method:: _create_validation_data_dict_for_genetic_algorithm() -> typing.Union[dict, pandas.DataFrame]
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._create_validation_data_dict_for_genetic_algorithm

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._create_validation_data_dict_for_genetic_algorithm

   .. py:method:: _get_validation_data_to_validate(valid_data: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData) -> pandas.DataFrame
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._get_validation_data_to_validate

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._get_validation_data_to_validate

   .. py:method:: _init_validation_df(bin_information: list = None, substrate_uptake_ids: list = None) -> dict
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._init_validation_df

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._init_validation_df

   .. py:method:: _correct_upper_and_lower_ranges_validation_data_for_sign(validation_range: list)
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._correct_upper_and_lower_ranges_validation_data_for_sign

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._correct_upper_and_lower_ranges_validation_data_for_sign

   .. py:method:: _init_results_objects()
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._init_results_objects

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._init_results_objects

   .. py:method:: _select_topn_enzymes(esc_results_df: pandas.DataFrame, nmbr_kcats_to_pick: int) -> pandas.DataFrame
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._select_topn_enzymes

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._select_topn_enzymes

   .. py:method:: _calculate_esc_variability(esc_df: pandas.DataFrame, enzyme_id: str) -> float
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._calculate_esc_variability

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._calculate_esc_variability

   .. py:method:: _esc_variability_larger_than_threshold(esc_variability: float) -> bool
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._esc_variability_larger_than_threshold

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._esc_variability_larger_than_threshold

   .. py:method:: _get_random_enzymes_to_evaluate()
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._get_random_enzymes_to_evaluate

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._get_random_enzymes_to_evaluate

   .. py:method:: _calculate_error_for_reactions(substrate_uptake_id: str, validation_df: pandas.DataFrame, reactions_to_validate: list, bin_id: typing.Union[float, int] = None) -> float
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._calculate_error_for_reactions

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._calculate_error_for_reactions

   .. py:method:: _parse_enzymes_to_evaluate(esc_topn_df: pandas.DataFrame) -> PAMparametrizer.PAM_parametrizer.pam_parametrizer.Reaction2KcatDict
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._parse_enzymes_to_evaluate

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._parse_enzymes_to_evaluate

   .. py:method:: _get_genetic_algorithm_json_files(subset: str = '')
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._get_genetic_algorithm_json_files

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._get_genetic_algorithm_json_files

   .. py:method:: _determine_enzymes_to_evaluate_for_all_bins(nmbr_kcats_to_pick) -> dict
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._determine_enzymes_to_evaluate_for_all_bins

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._determine_enzymes_to_evaluate_for_all_bins

   .. py:method:: _get_mutated_kcat_values_from_genetic_algorithm(results_filename: str = None)
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._get_mutated_kcat_values_from_genetic_algorithm

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._get_mutated_kcat_values_from_genetic_algorithm

   .. py:method:: _change_kcat_value_for_enzyme(enzyme_id: str, kcat_dict: dict) -> None
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._change_kcat_value_for_enzyme

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._change_kcat_value_for_enzyme

   .. py:method:: _remove_result_files(file_base: typing.Union[list, str]) -> None
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._remove_result_files

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._remove_result_files

   .. py:method:: _set_total_protein_constraint_to_equality()
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._set_total_protein_constraint_to_equality

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._set_total_protein_constraint_to_equality

   .. py:method:: _parse_row_with_enz_rxn_kcat_for_saving(enz_rxn_kcat_row: pandas.Series) -> list
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._parse_row_with_enz_rxn_kcat_for_saving

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._parse_row_with_enz_rxn_kcat_for_saving

   .. py:method:: _error_is_converging()
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._error_is_converging

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer._error_is_converging

   .. py:method:: plot_valid_data()
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.plot_valid_data

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.plot_valid_data

   .. py:method:: plot_simulation(fig, axs, return_fluxes: bool = False, save_esc=False, color: int = None, cbar_label: str = 'Iteration', sensitivity=True) -> matplotlib.pyplot.Figure
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.plot_simulation

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.plot_simulation

   .. py:method:: run_simulations_to_plot(substrate_uptake_id: str, substrate_rates: typing.Union[numpy.array, list, pandas.Series] = None, save_fluxes_esc: bool = False, sensitivity: bool = True) -> typing.Tuple[list, dict[float, float]]
      :canonical: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.run_simulations_to_plot

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.pam_parametrizer.PAMParametrizer.run_simulations_to_plot
