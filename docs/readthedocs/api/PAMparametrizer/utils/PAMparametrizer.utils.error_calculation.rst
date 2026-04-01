:py:mod:`PAMparametrizer.utils.error_calculation`
=================================================

.. py:module:: PAMparametrizer.utils.error_calculation

.. autodoc2-docstring:: PAMparametrizer.utils.error_calculation
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`calculate_r_squared_for_reaction <PAMparametrizer.utils.error_calculation.calculate_r_squared_for_reaction>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.error_calculation.calculate_r_squared_for_reaction
          :summary:
   * - :py:obj:`calculate_smape_for_reaction <PAMparametrizer.utils.error_calculation.calculate_smape_for_reaction>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.error_calculation.calculate_smape_for_reaction
          :summary:
   * - :py:obj:`custom_error <PAMparametrizer.utils.error_calculation.custom_error>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.error_calculation.custom_error
          :summary:
   * - :py:obj:`nanaverage <PAMparametrizer.utils.error_calculation.nanaverage>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.error_calculation.nanaverage
          :summary:
   * - :py:obj:`calculate_symmetric_mean_absolute_percentage_error <PAMparametrizer.utils.error_calculation.calculate_symmetric_mean_absolute_percentage_error>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.error_calculation.calculate_symmetric_mean_absolute_percentage_error
          :summary:

API
~~~

.. py:function:: calculate_r_squared_for_reaction(reaction_id: str, validation_data: pandas.DataFrame, substrate_uptake_id: str, fluxes: pandas.DataFrame) -> float
   :canonical: PAMparametrizer.utils.error_calculation.calculate_r_squared_for_reaction

   .. autodoc2-docstring:: PAMparametrizer.utils.error_calculation.calculate_r_squared_for_reaction

.. py:function:: calculate_smape_for_reaction(reaction_id: str, validation_data: pandas.DataFrame, substrate_uptake_id: str, fluxes: pandas.DataFrame) -> float
   :canonical: PAMparametrizer.utils.error_calculation.calculate_smape_for_reaction

   .. autodoc2-docstring:: PAMparametrizer.utils.error_calculation.calculate_smape_for_reaction

.. py:function:: custom_error(observed, simulated, lambda_factor=1.0)
   :canonical: PAMparametrizer.utils.error_calculation.custom_error

   .. autodoc2-docstring:: PAMparametrizer.utils.error_calculation.custom_error

.. py:function:: nanaverage(data: typing.Union[list], weights: dict = None, axis: int = None) -> typing.Iterable
   :canonical: PAMparametrizer.utils.error_calculation.nanaverage

   .. autodoc2-docstring:: PAMparametrizer.utils.error_calculation.nanaverage

.. py:function:: calculate_symmetric_mean_absolute_percentage_error(y_true: typing.Iterable[float], y_pred: typing.Iterable[float]) -> float
   :canonical: PAMparametrizer.utils.error_calculation.calculate_symmetric_mean_absolute_percentage_error

   .. autodoc2-docstring:: PAMparametrizer.utils.error_calculation.calculate_symmetric_mean_absolute_percentage_error
