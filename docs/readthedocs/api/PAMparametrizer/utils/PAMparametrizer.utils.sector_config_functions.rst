:py:mod:`PAMparametrizer.utils.sector_config_functions`
=======================================================

.. py:module:: PAMparametrizer.utils.sector_config_functions

.. autodoc2-docstring:: PAMparametrizer.utils.sector_config_functions
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`perform_linear_regression <PAMparametrizer.utils.sector_config_functions.perform_linear_regression>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.sector_config_functions.perform_linear_regression
          :summary:
   * - :py:obj:`reset_translational_sector <PAMparametrizer.utils.sector_config_functions.reset_translational_sector>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.sector_config_functions.reset_translational_sector
          :summary:
   * - :py:obj:`get_model_simulations_vs_sector <PAMparametrizer.utils.sector_config_functions.get_model_simulations_vs_sector>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.sector_config_functions.get_model_simulations_vs_sector
          :summary:
   * - :py:obj:`run_simulations <PAMparametrizer.utils.sector_config_functions.run_simulations>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.sector_config_functions.run_simulations
          :summary:
   * - :py:obj:`plot_translational_protein_vs_mu <PAMparametrizer.utils.sector_config_functions.plot_translational_protein_vs_mu>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.sector_config_functions.plot_translational_protein_vs_mu
          :summary:
   * - :py:obj:`plot_unused_protein_vs_mu <PAMparametrizer.utils.sector_config_functions.plot_unused_protein_vs_mu>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.sector_config_functions.plot_unused_protein_vs_mu
          :summary:
   * - :py:obj:`change_sector_parameters_with_config_dict <PAMparametrizer.utils.sector_config_functions.change_sector_parameters_with_config_dict>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.sector_config_functions.change_sector_parameters_with_config_dict
          :summary:
   * - :py:obj:`get_protein_sector_config <PAMparametrizer.utils.sector_config_functions.get_protein_sector_config>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.sector_config_functions.get_protein_sector_config
          :summary:
   * - :py:obj:`change_proteinsector_relation_from_growth_to_substrate_uptake <PAMparametrizer.utils.sector_config_functions.change_proteinsector_relation_from_growth_to_substrate_uptake>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.sector_config_functions.change_proteinsector_relation_from_growth_to_substrate_uptake
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`SectorParameterDict <PAMparametrizer.utils.sector_config_functions.SectorParameterDict>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.sector_config_functions.SectorParameterDict
          :summary:

API
~~~

.. py:data:: SectorParameterDict
   :canonical: PAMparametrizer.utils.sector_config_functions.SectorParameterDict
   :value: None

   .. autodoc2-docstring:: PAMparametrizer.utils.sector_config_functions.SectorParameterDict

.. py:function:: perform_linear_regression(x: typing.Iterable[typing.Union[float, int]], y: typing.Iterable[typing.Union[float, int]]) -> typing.Tuple[float, float]
   :canonical: PAMparametrizer.utils.sector_config_functions.perform_linear_regression

   .. autodoc2-docstring:: PAMparametrizer.utils.sector_config_functions.perform_linear_regression

.. py:function:: reset_translational_sector(pamodel: PAModelpy.PAModel, slope: float, intercept: float, new_id: str = None) -> PAModelpy.PAModel
   :canonical: PAMparametrizer.utils.sector_config_functions.reset_translational_sector

   .. autodoc2-docstring:: PAMparametrizer.utils.sector_config_functions.reset_translational_sector

.. py:function:: get_model_simulations_vs_sector(pamodel: PAModelpy.PAModel, sub_uptake_rxn: str, rxn_id_to_relate_to: str, substrate_range: typing.Iterable, intercept: float, slope: float, sector_name: str = 'translational_protein', to_save: str = None) -> pandas.DataFrame
   :canonical: PAMparametrizer.utils.sector_config_functions.get_model_simulations_vs_sector

   .. autodoc2-docstring:: PAMparametrizer.utils.sector_config_functions.get_model_simulations_vs_sector

.. py:function:: run_simulations(pamodel: PAModelpy.PAModel, substrate_rates: typing.Iterable[float], sub_uptake_id='EX_glc__D_e') -> list
   :canonical: PAMparametrizer.utils.sector_config_functions.run_simulations

   .. autodoc2-docstring:: PAMparametrizer.utils.sector_config_functions.run_simulations

.. py:function:: plot_translational_protein_vs_mu(literature: pandas.DataFrame, results: pandas.DataFrame, protein_fraction: float, measured_protein_fraction: float, oxygen_results: pandas.DataFrame = None, oxygen_rxn_id: str = None, return_fig: bool = False, configuration: PAModelpy.Config = None, literature_label: str = 'Schmidt et al (2016)', model_label: str = 'new iML1515 PAM', fig: matplotlib.pyplot.Figure = None, ax: matplotlib.pyplot.Axes = None) -> None
   :canonical: PAMparametrizer.utils.sector_config_functions.plot_translational_protein_vs_mu

   .. autodoc2-docstring:: PAMparametrizer.utils.sector_config_functions.plot_translational_protein_vs_mu

.. py:function:: plot_unused_protein_vs_mu(results: pandas.DataFrame, biomass_rxn: str) -> None
   :canonical: PAMparametrizer.utils.sector_config_functions.plot_unused_protein_vs_mu

   .. autodoc2-docstring:: PAMparametrizer.utils.sector_config_functions.plot_unused_protein_vs_mu

.. py:function:: change_sector_parameters_with_config_dict(pamodel: PAModelpy.PAModel, sector_config: PAMparametrizer.utils.sector_config_functions.SectorParameterDict, substrate_uptake_id: str, sector_id: str = 'TranslationalProteinSector') -> PAModelpy.PAModel
   :canonical: PAMparametrizer.utils.sector_config_functions.change_sector_parameters_with_config_dict

   .. autodoc2-docstring:: PAMparametrizer.utils.sector_config_functions.change_sector_parameters_with_config_dict

.. py:function:: get_protein_sector_config(pamodel: PAModelpy.PAModel, substrate_id: str, substrate_range: typing.Iterable[typing.Union[int, float]], rxn_to_relate_to: str = None, protein_sector: str = 'TranslationalProteinSector') -> PAMparametrizer.utils.sector_config_functions.SectorParameterDict
   :canonical: PAMparametrizer.utils.sector_config_functions.get_protein_sector_config

   .. autodoc2-docstring:: PAMparametrizer.utils.sector_config_functions.get_protein_sector_config

.. py:function:: change_proteinsector_relation_from_growth_to_substrate_uptake(pamodel: PAModelpy.PAModel, params: PAMparametrizer.utils.sector_config_functions.SectorParameterDict, sector_id: str, substrate_uptake_id: str = 'EX_glc__D_e', substrate_range: typing.Iterable[typing.Union[int, float]] = np.arange(-4, 0, 1), sector_name: str = 'unused_enzymes') -> PAMparametrizer.utils.sector_config_functions.SectorParameterDict
   :canonical: PAMparametrizer.utils.sector_config_functions.change_proteinsector_relation_from_growth_to_substrate_uptake

   .. autodoc2-docstring:: PAMparametrizer.utils.sector_config_functions.change_proteinsector_relation_from_growth_to_substrate_uptake
