:py:mod:`PAMparametrizer.PAM_parametrizer.KcatConstraintConfig`
===============================================================

.. py:module:: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig

.. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig
   :allowtitles:

Module Contents
---------------

Classes
~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`KcatConstraintConfigTable <PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable>`
     - .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable
          :summary:

API
~~~

.. py:class:: KcatConstraintConfigTable(df: typing.Optional[pandas.DataFrame] = None)
   :canonical: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable

   .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable

   .. rubric:: Initialization

   .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable.__init__

   .. py:attribute:: REQUIRED_COLUMNS
      :canonical: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable.REQUIRED_COLUMNS
      :value: ['enzyme_id', 'reaction_id', 'direction', 'min_kcat', 'max_kcat']

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable.REQUIRED_COLUMNS

   .. py:attribute:: DIFFUSION_LIMIT
      :canonical: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable.DIFFUSION_LIMIT
      :value: 1000000.0

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable.DIFFUSION_LIMIT

   .. py:attribute:: MIN_KCAT
      :canonical: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable.MIN_KCAT
      :value: 1e-06

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable.MIN_KCAT

   .. py:property:: df_model_constraints
      :canonical: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable.df_model_constraints
      :type: pandas.DataFrame

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable.df_model_constraints

   .. py:method:: _validate_input_df(df: pandas.DataFrame)
      :canonical: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable._validate_input_df

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable._validate_input_df

   .. py:method:: get(enzyme_id: str, reaction_id: str, direction: str) -> dict
      :canonical: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable.get

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable.get

   .. py:method:: get_in_model_constraints(enzyme_id: str, reaction_id: str, direction: str) -> dict
      :canonical: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable.get_in_model_constraints

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable.get_in_model_constraints

   .. py:method:: add(enzyme_id: str, reaction_id: str, direction: str, min_kcat: typing.Optional[float] = MIN_KCAT, max_kcat: typing.Optional[float] = DIFFUSION_LIMIT) -> dict
      :canonical: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable.add

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable.add

   .. py:method:: has_constraint(enzyme_id: str, reaction_id: str, direction: str) -> bool
      :canonical: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable.has_constraint

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable.has_constraint

   .. py:method:: ensure_kcats_in_pam_info_file_are_within_bounds(pam_info_file: str) -> None
      :canonical: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable.ensure_kcats_in_pam_info_file_are_within_bounds

      .. autodoc2-docstring:: PAMparametrizer.PAM_parametrizer.KcatConstraintConfig.KcatConstraintConfigTable.ensure_kcats_in_pam_info_file_are_within_bounds
