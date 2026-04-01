:py:mod:`PAMparametrizer.utils.pamparametrizer_setup`
=====================================================

.. py:module:: PAMparametrizer.utils.pamparametrizer_setup

.. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_setup
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`save_sector_information_to_excel <PAMparametrizer.utils.pamparametrizer_setup.save_sector_information_to_excel>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_setup.save_sector_information_to_excel
          :summary:
   * - :py:obj:`_get_pam_parameter_information_from_excel <PAMparametrizer.utils.pamparametrizer_setup._get_pam_parameter_information_from_excel>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_setup._get_pam_parameter_information_from_excel
          :summary:
   * - :py:obj:`set_up_sector_config <PAMparametrizer.utils.pamparametrizer_setup.set_up_sector_config>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_setup.set_up_sector_config
          :summary:
   * - :py:obj:`set_up_sector_config_from_diagnostic_file <PAMparametrizer.utils.pamparametrizer_setup.set_up_sector_config_from_diagnostic_file>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_setup.set_up_sector_config_from_diagnostic_file
          :summary:
   * - :py:obj:`get_kcat_constraints <PAMparametrizer.utils.pamparametrizer_setup.get_kcat_constraints>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_setup.get_kcat_constraints
          :summary:

API
~~~

.. py:function:: save_sector_information_to_excel(param_vs_lin_rxn: PAMparametrizer.utils.sector_config_functions.SectorParameterDict, lin_rxn_id: str, sector_id: typing.Literal[Translational, UnusedEnzyme], biomass_rxn: str = None, param_vs_growth: PAMparametrizer.utils.sector_config_functions.SectorParameterDict = None, pam_data_file: str = None, output_file_path: str = None, substrate_range: typing.Iterable[typing.Union[float, int]] = np.arange(-4, 0, 1), model_name: str = 'iML1515') -> None
   :canonical: PAMparametrizer.utils.pamparametrizer_setup.save_sector_information_to_excel

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_setup.save_sector_information_to_excel

.. py:function:: _get_pam_parameter_information_from_excel(pam_data_file: str, model_name: str, sector_id: str, output_file_path: str = None) -> typing.Tuple[dict[str, pandas.DataFrame], str]
   :canonical: PAMparametrizer.utils.pamparametrizer_setup._get_pam_parameter_information_from_excel

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_setup._get_pam_parameter_information_from_excel

.. py:function:: set_up_sector_config(pam_info_file: str, sectors_not_related_to_growth: typing.List[str]) -> typing.Dict[str, PAMparametrizer.PAM_parametrizer.SectorConfig]
   :canonical: PAMparametrizer.utils.pamparametrizer_setup.set_up_sector_config

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_setup.set_up_sector_config

.. py:function:: set_up_sector_config_from_diagnostic_file(diagnostic_file: str, substrate_uptake_id: str = 'EX_glc__D_e') -> typing.Dict[str, typing.Dict[typing.Literal[slope, intercept], float]]
   :canonical: PAMparametrizer.utils.pamparametrizer_setup.set_up_sector_config_from_diagnostic_file

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_setup.set_up_sector_config_from_diagnostic_file

.. py:function:: get_kcat_constraints(pam_info_file: str) -> PAMparametrizer.PAM_parametrizer.KcatConstraintConfigTable
   :canonical: PAMparametrizer.utils.pamparametrizer_setup.get_kcat_constraints

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_setup.get_kcat_constraints
