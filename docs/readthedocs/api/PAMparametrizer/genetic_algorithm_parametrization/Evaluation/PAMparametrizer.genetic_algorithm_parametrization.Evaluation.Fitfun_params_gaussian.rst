:py:mod:`PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian`
=============================================================================================

.. py:module:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian

.. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian
   :allowtitles:

Module Contents
---------------

Classes
~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`FitnessEvaluation <PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation>`
     - .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`FILE_PATH <PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FILE_PATH>`
     - .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FILE_PATH
          :summary:
   * - :py:obj:`DATA_PATH <PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.DATA_PATH>`
     - .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.DATA_PATH
          :summary:
   * - :py:obj:`DIFUSSIONLIMIT <PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.DIFUSSIONLIMIT>`
     - .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.DIFUSSIONLIMIT
          :summary:

API
~~~

.. py:data:: FILE_PATH
   :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FILE_PATH
   :value: 'Path(...)'

   .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FILE_PATH

.. py:data:: DATA_PATH
   :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.DATA_PATH
   :value: 'joinpath(...)'

   .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.DATA_PATH

.. py:data:: DIFUSSIONLIMIT
   :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.DIFUSSIONLIMIT
   :value: None

   .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.DIFUSSIONLIMIT

.. py:class:: FitnessEvaluation(model=None, sector_configs_per_substrate: dict = None, fixed_attr_list=[], objective_id=str(), valid_data=dict(), sigma_denominator: int = 10, error_weights: dict = {}, substrate_uptake_rates={'EX_glc__D_e': [0.7, 11.3]}, substrate_uptake_id='EX_glc__D_e')
   :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation

   .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation

   .. rubric:: Initialization

   .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation.__init__

   .. py:attribute:: KCAT_MU
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation.KCAT_MU
      :value: 13.7

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation.KCAT_MU

   .. py:attribute:: KCAT_SIGMA
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation.KCAT_SIGMA
      :value: 100.0

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation.KCAT_SIGMA

   .. py:attribute:: NUM_KCATS
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation.NUM_KCATS
      :value: 5

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation.NUM_KCATS

   .. py:method:: init_model(model=None)
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation.init_model

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation.init_model

   .. py:method:: compute_individuals_properties(pop) -> list
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation.compute_individuals_properties

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation.compute_individuals_properties

   .. py:method:: _determine_changed_kcats(individual)
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation._determine_changed_kcats

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation._determine_changed_kcats

   .. py:method:: attribute_generator(probability=0, kcat_list=[]) -> int
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation.attribute_generator

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation.attribute_generator

   .. py:method:: init_fitness() -> dict
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation.init_fitness

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation.init_fitness

   .. py:method:: init_attribute(enzymes_to_eval: list, direction: list, rxns_to_eval: list) -> dict
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation.init_attribute

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation.init_attribute

   .. py:method:: init_individual() -> dict
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation.init_individual

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation.init_individual

   .. py:method:: eval_fitness(individual) -> int
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation.eval_fitness

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation.eval_fitness

   .. py:method:: _reset_sectors()
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation._reset_sectors

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation._reset_sectors

   .. py:method:: _get_old_kcats(model: PAModelpy.PAModel, individual) -> typing.List[float]
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation._get_old_kcats

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation._get_old_kcats

   .. py:method:: _mutate_kcat_value(kcat: float, min_kcat: float = 1e-06, max_kcat: float = DIFUSSIONLIMIT, sensitivity: float = 0.5, toolbox: deap.base.Toolbox = None) -> float
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation._mutate_kcat_value

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation._mutate_kcat_value

   .. py:method:: _change_kcat_values_for_individual(individual, kcat_values: list = [])
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation._change_kcat_values_for_individual

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation._change_kcat_values_for_individual

   .. py:method:: _calculate_simulation_error(flux_df: pandas.DataFrame, substrate_reaction: str)
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation._calculate_simulation_error

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_gaussian.FitnessEvaluation._calculate_simulation_error
