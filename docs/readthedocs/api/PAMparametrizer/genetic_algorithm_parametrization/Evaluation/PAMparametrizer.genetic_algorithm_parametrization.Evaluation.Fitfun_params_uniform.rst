:py:mod:`PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_uniform`
============================================================================================

.. py:module:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_uniform

.. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_uniform
   :allowtitles:

Module Contents
---------------

Classes
~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`FitnessEvaluation <PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_uniform.FitnessEvaluation>`
     - .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_uniform.FitnessEvaluation
          :summary:

API
~~~

.. py:class:: FitnessEvaluation(*args, **kwargs)
   :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_uniform.FitnessEvaluation

   Bases: :py:obj:`PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_uniform.FitnessEvaluation`

   .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_uniform.FitnessEvaluation

   .. rubric:: Initialization

   .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_uniform.FitnessEvaluation.__init__

   .. py:method:: attribute_generator(probability: float = 0, kcat_list: typing.Optional[typing.List[float]] = [], kcat_bounds: typing.Optional[typing.List[dict]] = {}) -> int
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_uniform.FitnessEvaluation.attribute_generator

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_uniform.FitnessEvaluation.attribute_generator

   .. py:method:: _mutate_kcat_value(kcat: float, min_kcat: float = 0.001, max_kcat: float = DIFUSSIONLIMIT, sensitivity: float = 0.5, toolbox: deap.base.Toolbox = None) -> float
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_uniform.FitnessEvaluation._mutate_kcat_value

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_uniform.FitnessEvaluation._mutate_kcat_value

   .. py:method:: _mut_kcat_uniform(kcat_list: typing.List[float], indpb: float, min_kcat: float = 1e-06, max_kcat: float = DIFUSSIONLIMIT) -> typing.List[float]
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_uniform.FitnessEvaluation._mut_kcat_uniform

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.Evaluation.Fitfun_params_uniform.FitnessEvaluation._mut_kcat_uniform
