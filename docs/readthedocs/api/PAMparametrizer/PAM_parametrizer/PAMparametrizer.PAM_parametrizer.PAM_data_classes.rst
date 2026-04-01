:py:mod:`PAMparametrizer.PAM_parametrizer.PAM_data_classes`
===========================================================

.. py:module:: PAMparametrizer.PAM_parametrizer.PAM_data_classes

.. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes
   :allowtitles:

Module Contents
---------------

Classes
~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`SectorConfig <PAMparametrizer.PAM_parametrizer.PAM_data_classes.SectorConfig>`
     - .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.SectorConfig
          :summary:
   * - :py:obj:`ValidationData <PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData>`
     - .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData
          :summary:
   * - :py:obj:`HyperParameters <PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters>`
     - .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters
          :summary:
   * - :py:obj:`ParametrizationResults <PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults>`
     - .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults
          :summary:
   * - :py:obj:`FluxResults <PAMparametrizer.PAM_parametrizer.PAM_data_classes.FluxResults>`
     - .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.FluxResults
          :summary:

API
~~~

.. py:class:: SectorConfig()
   :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.SectorConfig

   Bases: :py:obj:`typing.TypedDict`

   .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.SectorConfig

   .. rubric:: Initialization

   .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.SectorConfig.__init__

   .. py:attribute:: sectorname
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.SectorConfig.sectorname
      :type: str
      :value: None

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.SectorConfig.sectorname

   .. py:attribute:: slope
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.SectorConfig.slope
      :type: float
      :value: None

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.SectorConfig.slope

   .. py:attribute:: intercept
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.SectorConfig.intercept
      :type: float
      :value: None

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.SectorConfig.intercept

   .. py:attribute:: substrate_range
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.SectorConfig.substrate_range
      :type: typing.Iterable[typing.Union[float, int]]
      :value: None

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.SectorConfig.substrate_range

.. py:class:: ValidationData
   :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData

   .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData

   .. py:attribute:: valid_data
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData.valid_data
      :type: pandas.DataFrame
      :value: None

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData.valid_data

   .. py:attribute:: id
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData.id
      :type: str
      :value: None

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData.id

   .. py:attribute:: validation_range
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData.validation_range
      :type: list
      :value: None

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData.validation_range

   .. py:attribute:: sampled_valid_data
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData.sampled_valid_data
      :type: pandas.DataFrame
      :value: None

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData.sampled_valid_data

   .. py:attribute:: inactive_exchanges
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData.inactive_exchanges
      :type: typing.Union[typing.List, None]
      :value: 'field(...)'

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData.inactive_exchanges

   .. py:attribute:: sector_configs
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData.sector_configs
      :type: dict
      :value: 'field(...)'

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData.sector_configs

   .. py:attribute:: _reactions_to_validate
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData._reactions_to_validate
      :type: list
      :value: 'field(...)'

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData._reactions_to_validate

   .. py:attribute:: biomass_reaction_extension
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData.biomass_reaction_extension
      :type: str
      :value: 'BIOMASS'

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData.biomass_reaction_extension

   .. py:attribute:: exchange_reaction_extension
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData.exchange_reaction_extension
      :type: str
      :value: 'EX'

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData.exchange_reaction_extension

   .. py:attribute:: _reactions_to_plot
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData._reactions_to_plot
      :value: None

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData._reactions_to_plot

   .. py:method:: __post_init__()
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData.__post_init__

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData.__post_init__

   .. py:method:: _add_inactive_exchanges_to_validation_df() -> None
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData._add_inactive_exchanges_to_validation_df

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData._add_inactive_exchanges_to_validation_df

   .. py:method:: _get_biomass_reactions() -> list
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData._get_biomass_reactions

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData._get_biomass_reactions

   .. py:property:: biomass_reactions
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData.biomass_reactions
      :type: list

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData.biomass_reactions

   .. py:method:: _get_reactions_to_validate() -> list
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData._get_reactions_to_validate

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData._get_reactions_to_validate

   .. py:property:: reactions_to_validate
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData.reactions_to_validate
      :type: list

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ValidationData.reactions_to_validate

.. py:class:: HyperParameters
   :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters

   .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters

   .. py:attribute:: threshold_error
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.threshold_error
      :value: 0.9

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.threshold_error

   .. py:attribute:: threshold_iteration
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.threshold_iteration
      :value: 100

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.threshold_iteration

   .. py:attribute:: threshold_convergence
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.threshold_convergence
      :value: 0.09

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.threshold_convergence

   .. py:attribute:: threshold_nmbr_convergence
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.threshold_nmbr_convergence
      :value: 2

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.threshold_nmbr_convergence

   .. py:attribute:: number_of_kcats_to_mutate
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.number_of_kcats_to_mutate
      :value: 5

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.number_of_kcats_to_mutate

   .. py:attribute:: number_of_bins
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.number_of_bins
      :value: 5

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.number_of_bins

   .. py:attribute:: bin_resolution
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.bin_resolution
      :value: 5

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.bin_resolution

   .. py:attribute:: bin_split_deviation_threshold
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.bin_split_deviation_threshold
      :value: 0.2

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.bin_split_deviation_threshold

   .. py:attribute:: genetic_algorithm_filename_base
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.genetic_algorithm_filename_base
      :value: 'genetic_algorithm_run_'

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.genetic_algorithm_filename_base

   .. py:attribute:: filename_extension
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.filename_extension
      :value: <Multiline-String>

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.filename_extension

   .. py:attribute:: genetic_algorithm
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.genetic_algorithm
      :type: typing.Callable
      :value: None

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.genetic_algorithm

   .. py:attribute:: genetic_algorithm_hyperparams
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.genetic_algorithm_hyperparams
      :value: None

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.HyperParameters.genetic_algorithm_hyperparams

.. py:class:: ParametrizationResults
   :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults

   .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults

   .. py:attribute:: substrate_uptake_reactions
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.substrate_uptake_reactions
      :type: list
      :value: None

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.substrate_uptake_reactions

   .. py:attribute:: esc_df
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.esc_df
      :type: pandas.DataFrame
      :value: 'field(...)'

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.esc_df

   .. py:attribute:: sensitive_enzymes
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.sensitive_enzymes
      :type: pandas.DataFrame
      :value: 'field(...)'

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.sensitive_enzymes

   .. py:attribute:: flux_results
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.flux_results
      :type: cobra.DictList
      :value: None

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.flux_results

   .. py:attribute:: _color
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults._color
      :type: int
      :value: 440154

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults._color

   .. py:attribute:: best_individuals
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.best_individuals
      :type: pandas.DataFrame
      :value: 'field(...)'

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.best_individuals

   .. py:attribute:: computational_time
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.computational_time
      :type: pandas.DataFrame
      :value: 'field(...)'

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.computational_time

   .. py:attribute:: final_errors
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.final_errors
      :type: pandas.DataFrame
      :value: 'field(...)'

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.final_errors

   .. py:method:: initiate_flux_results()
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.initiate_flux_results

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.initiate_flux_results

   .. py:method:: initiate_result_dfs(reactions_to_validate: dict) -> None
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.initiate_result_dfs

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.initiate_result_dfs

   .. py:method:: initiate_bins_to_change() -> None
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.initiate_bins_to_change

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.initiate_bins_to_change

   .. py:method:: add_new_substrate_source(substrate_uptake_id: str, reactions_to_validate: list) -> None
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.add_new_substrate_source

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.add_new_substrate_source

   .. py:method:: add_fluxes(pamodel, bin_id: typing.Union[float, int, str], substrate_reaction_id: str, substrate_uptake_rate: typing.Union[float, int], fluxes_abs: bool = True)
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.add_fluxes

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.add_fluxes

   .. py:method:: add_fluxes_from_fluxdict(flux_dict: dict, bin_id: typing.Union[float, int, str], substrate_reaction_id: str, substrate_uptake_rate: typing.Union[float, int], fluxes_abs: bool = True)
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.add_fluxes_from_fluxdict

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.add_fluxes_from_fluxdict

   .. py:method:: add_enzyme_sensitivity_coefficients(enzyme_sensitivity_coeff: pandas.DataFrame, bin_id: typing.Union[float, int, str], substrate_uptake_rate: typing.Union[float, int])
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.add_enzyme_sensitivity_coefficients

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.add_enzyme_sensitivity_coefficients

   .. py:method:: add_error_to_error_df(substrate_uptake_id: str, bin_id: typing.Union[int, str], error: float)
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.add_error_to_error_df

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.add_error_to_error_df

   .. py:method:: add_best_individuals(run_id: str, best_indiv_enz_rxn_kcat: list, ga_error: float) -> None
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.add_best_individuals

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.add_best_individuals

   .. py:method:: add_computational_time(run_id: str, time_in_sec: float)
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.add_computational_time

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.add_computational_time

   .. py:method:: add_final_error(run_id: str, final_error: float)
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.add_final_error

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.add_final_error

   .. py:method:: remove_simulations_from_flux_df(substrate_uptake_id, bin_id)
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.remove_simulations_from_flux_df

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.ParametrizationResults.remove_simulations_from_flux_df

.. py:class:: FluxResults
   :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.FluxResults

   .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.FluxResults

   .. py:attribute:: id
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.FluxResults.id
      :type: str
      :value: None

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.FluxResults.id

   .. py:attribute:: error_df
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.FluxResults.error_df
      :type: pandas.DataFrame
      :value: None

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.FluxResults.error_df

   .. py:attribute:: fluxes_df
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.FluxResults.fluxes_df
      :type: pandas.DataFrame
      :value: None

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.FluxResults.fluxes_df

   .. py:attribute:: substrate_range
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.FluxResults.substrate_range
      :type: list
      :value: 'field(...)'

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.FluxResults.substrate_range

   .. py:method:: initiate_result_dfs(reactions_to_validate: list) -> None
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.FluxResults.initiate_result_dfs

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.FluxResults.initiate_result_dfs

   .. py:method:: add_fluxes(pamodel, bin_id: typing.Union[float, int, str], substrate_uptake_rate: typing.Union[float, int], fluxes_abs: bool = True)
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.FluxResults.add_fluxes

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.FluxResults.add_fluxes

   .. py:method:: add_fluxes_from_fluxdict(flux_dict: dict, bin_id: typing.Union[float, int, str], substrate_uptake_rate: typing.Union[float, int], fluxes_abs: bool = True)
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.FluxResults.add_fluxes_from_fluxdict

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.FluxResults.add_fluxes_from_fluxdict

   .. py:method:: add_error_to_error_df(bin_id: typing.Union[int, str], error: float)
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.FluxResults.add_error_to_error_df

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.FluxResults.add_error_to_error_df

   .. py:method:: remove_simulations_from_flux_df(bin_id) -> None
      :canonical: PAMparametrizer.PAM_parametrizer.PAM_data_classes.FluxResults.remove_simulations_from_flux_df

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.PAM_data_classes.FluxResults.remove_simulations_from_flux_df
