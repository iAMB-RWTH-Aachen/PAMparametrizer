:py:mod:`PAMparametrizer.utils.preprocessing`
=============================================

.. py:module:: PAMparametrizer.utils.preprocessing

.. autodoc2-docstring:: PAMparametrizer.utils.preprocessing
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`extract_locus_tags <PAMparametrizer.utils.preprocessing.extract_locus_tags>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.extract_locus_tags
          :summary:
   * - :py:obj:`create_id_mapper_from_model <PAMparametrizer.utils.preprocessing.create_id_mapper_from_model>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.create_id_mapper_from_model
          :summary:
   * - :py:obj:`has_same_reactants_and_products <PAMparametrizer.utils.preprocessing.has_same_reactants_and_products>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.has_same_reactants_and_products
          :summary:
   * - :py:obj:`create_genetokeggid_mapper <PAMparametrizer.utils.preprocessing.create_genetokeggid_mapper>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.create_genetokeggid_mapper
          :summary:
   * - :py:obj:`replace_locustags_in_text <PAMparametrizer.utils.preprocessing.replace_locustags_in_text>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.replace_locustags_in_text
          :summary:
   * - :py:obj:`map_kcat_values_to_reaction_protein_association <PAMparametrizer.utils.preprocessing.map_kcat_values_to_reaction_protein_association>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.map_kcat_values_to_reaction_protein_association
          :summary:
   * - :py:obj:`assign_missing_gprs <PAMparametrizer.utils.preprocessing.assign_missing_gprs>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.assign_missing_gprs
          :summary:
   * - :py:obj:`assign_directionalities_for_kcat_relations <PAMparametrizer.utils.preprocessing.assign_directionalities_for_kcat_relations>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.assign_directionalities_for_kcat_relations
          :summary:
   * - :py:obj:`assign_defaults_for_proteins_without_mapping <PAMparametrizer.utils.preprocessing.assign_defaults_for_proteins_without_mapping>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.assign_defaults_for_proteins_without_mapping
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`DEFAULT_KCAT <PAMparametrizer.utils.preprocessing.DEFAULT_KCAT>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.DEFAULT_KCAT
          :summary:
   * - :py:obj:`DEFAULT_MOLMASS <PAMparametrizer.utils.preprocessing.DEFAULT_MOLMASS>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.DEFAULT_MOLMASS
          :summary:
   * - :py:obj:`DEFAULT_PROT_LENGTH <PAMparametrizer.utils.preprocessing.DEFAULT_PROT_LENGTH>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.DEFAULT_PROT_LENGTH
          :summary:

API
~~~

.. py:data:: DEFAULT_KCAT
   :canonical: PAMparametrizer.utils.preprocessing.DEFAULT_KCAT
   :value: 13.7

   .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.DEFAULT_KCAT

.. py:data:: DEFAULT_MOLMASS
   :canonical: PAMparametrizer.utils.preprocessing.DEFAULT_MOLMASS
   :value: 39959.4825

   .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.DEFAULT_MOLMASS

.. py:data:: DEFAULT_PROT_LENGTH
   :canonical: PAMparametrizer.utils.preprocessing.DEFAULT_PROT_LENGTH
   :value: 325

   .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.DEFAULT_PROT_LENGTH

.. py:function:: extract_locus_tags(text: str, locustag_regex: str) -> typing.List[str]
   :canonical: PAMparametrizer.utils.preprocessing.extract_locus_tags

   .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.extract_locus_tags

.. py:function:: create_id_mapper_from_model(model: cobra.Model, rxn_annotation_keys: list[str] = ['kegg.reaction', 'ec-code'], met_annotation_key: str = 'kegg.compound', exclude: typing.List[str] = ['Growth', 'ATPM', 'BIOMASS']) -> pandas.DataFrame
   :canonical: PAMparametrizer.utils.preprocessing.create_id_mapper_from_model

   .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.create_id_mapper_from_model

.. py:function:: has_same_reactants_and_products(reaction)
   :canonical: PAMparametrizer.utils.preprocessing.has_same_reactants_and_products

   .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.has_same_reactants_and_products

.. py:function:: create_genetokeggid_mapper(model: cobra.Model) -> pandas.DataFrame
   :canonical: PAMparametrizer.utils.preprocessing.create_genetokeggid_mapper

   .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.create_genetokeggid_mapper

.. py:function:: replace_locustags_in_text(text: str, id_map: typing.Dict[str, str]) -> str
   :canonical: PAMparametrizer.utils.preprocessing.replace_locustags_in_text

   .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.replace_locustags_in_text

.. py:function:: map_kcat_values_to_reaction_protein_association(id_mapper: pandas.DataFrame, gotenzymes_df: pandas.DataFrame) -> pandas.DataFrame
   :canonical: PAMparametrizer.utils.preprocessing.map_kcat_values_to_reaction_protein_association

   .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.map_kcat_values_to_reaction_protein_association

.. py:function:: assign_missing_gprs(df: pandas.DataFrame, use_ec: bool = False) -> pandas.DataFrame
   :canonical: PAMparametrizer.utils.preprocessing.assign_missing_gprs

   .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.assign_missing_gprs

.. py:function:: assign_directionalities_for_kcat_relations(eco_enzymes: pandas.DataFrame) -> pandas.DataFrame
   :canonical: PAMparametrizer.utils.preprocessing.assign_directionalities_for_kcat_relations

   .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.assign_directionalities_for_kcat_relations

.. py:function:: assign_defaults_for_proteins_without_mapping(eco_enzymes: pandas.DataFrame, default_kcat: float = DEFAULT_KCAT, default_molmass: typing.Union[float, int] = DEFAULT_MOLMASS, default_protein_length: int = DEFAULT_PROT_LENGTH) -> pandas.DataFrame
   :canonical: PAMparametrizer.utils.preprocessing.assign_defaults_for_proteins_without_mapping

   .. autodoc2-docstring:: PAMparametrizer.utils.preprocessing.assign_defaults_for_proteins_without_mapping
