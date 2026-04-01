:py:mod:`PAMparametrizer.utils.pam_generation`
==============================================

.. py:module:: PAMparametrizer.utils.pam_generation

.. autodoc2-docstring:: PAMparametrizer.utils.pam_generation
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`create_pamodel_from_diagnostics_file <PAMparametrizer.utils.pam_generation.create_pamodel_from_diagnostics_file>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation.create_pamodel_from_diagnostics_file
          :summary:
   * - :py:obj:`get_rxn2kcat_protein2gene_dict <PAMparametrizer.utils.pam_generation.get_rxn2kcat_protein2gene_dict>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation.get_rxn2kcat_protein2gene_dict
          :summary:
   * - :py:obj:`_extract_reaction_id_from_catalytic_reaction_id <PAMparametrizer.utils.pam_generation._extract_reaction_id_from_catalytic_reaction_id>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation._extract_reaction_id_from_catalytic_reaction_id
          :summary:
   * - :py:obj:`_get_rxn2kcat_as_series <PAMparametrizer.utils.pam_generation._get_rxn2kcat_as_series>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation._get_rxn2kcat_as_series
          :summary:
   * - :py:obj:`search_index_in_parameter_file <PAMparametrizer.utils.pam_generation.search_index_in_parameter_file>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation.search_index_in_parameter_file
          :summary:
   * - :py:obj:`create_new_aes_parameter_file <PAMparametrizer.utils.pam_generation.create_new_aes_parameter_file>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation.create_new_aes_parameter_file
          :summary:
   * - :py:obj:`setup_yeast_pam <PAMparametrizer.utils.pam_generation.setup_yeast_pam>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation.setup_yeast_pam
          :summary:
   * - :py:obj:`set_up_yeast_config <PAMparametrizer.utils.pam_generation.set_up_yeast_config>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation.set_up_yeast_config
          :summary:
   * - :py:obj:`setup_pputida_pam <PAMparametrizer.utils.pam_generation.setup_pputida_pam>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation.setup_pputida_pam
          :summary:
   * - :py:obj:`setup_cglutamicum_pam <PAMparametrizer.utils.pam_generation.setup_cglutamicum_pam>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation.setup_cglutamicum_pam
          :summary:
   * - :py:obj:`turn_of_exchanges <PAMparametrizer.utils.pam_generation.turn_of_exchanges>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation.turn_of_exchanges
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`DEFAULT_MOLMASS <PAMparametrizer.utils.pam_generation.DEFAULT_MOLMASS>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation.DEFAULT_MOLMASS
          :summary:
   * - :py:obj:`DEFAULT_KCAT <PAMparametrizer.utils.pam_generation.DEFAULT_KCAT>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation.DEFAULT_KCAT
          :summary:

API
~~~

.. py:data:: DEFAULT_MOLMASS
   :canonical: PAMparametrizer.utils.pam_generation.DEFAULT_MOLMASS
   :value: 39959.4825

   .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation.DEFAULT_MOLMASS

.. py:data:: DEFAULT_KCAT
   :canonical: PAMparametrizer.utils.pam_generation.DEFAULT_KCAT
   :value: 11

   .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation.DEFAULT_KCAT

.. py:function:: create_pamodel_from_diagnostics_file(file_path: str, model: PAModelpy.PAModel, sheet_name: str = 'Best_Individuals', enzyme_sector_update: bool = True, other_enzyme_id_pattern: str = 'E[0-9][0-9]*|Enzyme*', substrate_uptake_id: str = 'EX_glc__D_e') -> PAModelpy.PAModel
   :canonical: PAMparametrizer.utils.pam_generation.create_pamodel_from_diagnostics_file

   .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation.create_pamodel_from_diagnostics_file

.. py:function:: get_rxn2kcat_protein2gene_dict(param_file_path: str, model_file_path: str) -> typing.Tuple[dict[str, dict[str, dict[typing.Literal[f, b, molmass, protein_reaction_association], float]]], dict[str, str]]
   :canonical: PAMparametrizer.utils.pam_generation.get_rxn2kcat_protein2gene_dict

   .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation.get_rxn2kcat_protein2gene_dict

.. py:function:: _extract_reaction_id_from_catalytic_reaction_id(input_str: str, default_enzyme_id_pattern: str = 'E[0-9][0-9]*|Enzyme_*') -> str
   :canonical: PAMparametrizer.utils.pam_generation._extract_reaction_id_from_catalytic_reaction_id

   .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation._extract_reaction_id_from_catalytic_reaction_id

.. py:function:: _get_rxn2kcat_as_series(rxn2kcat: dict[str, dict], name: str)
   :canonical: PAMparametrizer.utils.pam_generation._get_rxn2kcat_as_series

   .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation._get_rxn2kcat_as_series

.. py:function:: search_index_in_parameter_file(df: pandas.DataFrame, protein: str, reaction: str, direction: str) -> pandas.Index
   :canonical: PAMparametrizer.utils.pam_generation.search_index_in_parameter_file

   .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation.search_index_in_parameter_file

.. py:function:: create_new_aes_parameter_file(old_param_file: str, diagnostics_file_path: str, new_aes_file: str, diagnostics_sheet_name: str = 'Best_Individuals') -> None
   :canonical: PAMparametrizer.utils.pam_generation.create_new_aes_parameter_file

   .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation.create_new_aes_parameter_file

.. py:function:: setup_yeast_pam(pam_info_file: str = os.path.join('Results', '1_preprocessing', 'proteinAllocationModel_yeast9_EnzymaticData_TurnUp.xlsx'), model: str = 'Models/yeast9.xml', config: PAModelpy.configuration.Config = None, total_protein: typing.Union[bool, float] = 0.28697423725932236, active_enzymes: bool = True, translational_enzymes: bool = True, unused_enzymes: bool = True, sensitivity=True) -> PAModelpy.PAModel
   :canonical: PAMparametrizer.utils.pam_generation.setup_yeast_pam

   .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation.setup_yeast_pam

.. py:function:: set_up_yeast_config()
   :canonical: PAMparametrizer.utils.pam_generation.set_up_yeast_config

   .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation.set_up_yeast_config

.. py:function:: setup_pputida_pam(pam_info_file: str = os.path.join('Results', '1_preprocessing', 'proteinAllocationModel_iJN1463_EnzymaticData_250915.xlsx'), model: str = 'Models/iJN1463.xml', total_protein: typing.Union[bool, float] = 0.3, active_enzymes: bool = True, translational_enzymes: bool = True, unused_enzymes: bool = True, sensitivity=True)
   :canonical: PAMparametrizer.utils.pam_generation.setup_pputida_pam

   .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation.setup_pputida_pam

.. py:function:: setup_cglutamicum_pam(pam_info_file: str = os.path.join('Results', '1_preprocessing', 'proteinAllocationModel_iCGB21FR_EnzymaticData_250915.xlsx'), model: str = 'Models/iCGB21FR_annotated_copyable.xml', total_protein: typing.Union[bool, float] = 0.3, active_enzymes: bool = True, translational_enzymes: bool = True, unused_enzymes: bool = True, sensitivity=True)
   :canonical: PAMparametrizer.utils.pam_generation.setup_cglutamicum_pam

   .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation.setup_cglutamicum_pam

.. py:function:: turn_of_exchanges(model: PAModelpy.PAModel, exchanges_to_exclude: typing.List[str]) -> None
   :canonical: PAMparametrizer.utils.pam_generation.turn_of_exchanges

   .. autodoc2-docstring:: PAMparametrizer.utils.pam_generation.turn_of_exchanges
