:py:mod:`PAMparametrizer.genetic_algorithm_parametrization.ga_param`
====================================================================

.. py:module:: PAMparametrizer.genetic_algorithm_parametrization.ga_param

.. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.ga_param
   :allowtitles:

Module Contents
---------------

Classes
~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`Genetic_Algorithm <PAMparametrizer.genetic_algorithm_parametrization.ga_param.Genetic_Algorithm>`
     - .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.ga_param.Genetic_Algorithm
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`print_time <PAMparametrizer.genetic_algorithm_parametrization.ga_param.print_time>`
     - .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.ga_param.print_time
          :summary:

API
~~~

.. py:data:: print_time
   :canonical: PAMparametrizer.genetic_algorithm_parametrization.ga_param.print_time
   :value: None

   .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.ga_param.print_time

.. py:class:: Genetic_Algorithm(crossover_probability: float = 0.8, mutation_probability: float = 0.5, number_generations: int = 20, time_limit: int = 600)
   :canonical: PAMparametrizer.genetic_algorithm_parametrization.ga_param.Genetic_Algorithm

   .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.ga_param.Genetic_Algorithm

   .. rubric:: Initialization

   .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.ga_param.Genetic_Algorithm.__init__

   .. py:method:: init_pop(toolbox, population_size, evaluate_fitness=True) -> list
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.ga_param.Genetic_Algorithm.init_pop

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.ga_param.Genetic_Algorithm.init_pop

   .. py:method:: evaluate_pop(pop: list, toolbox: deap.base.Toolbox) -> list
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.ga_param.Genetic_Algorithm.evaluate_pop

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.ga_param.Genetic_Algorithm.evaluate_pop

   .. py:method:: main(pop: list, toolbox: deap.base.Toolbox, start_time: float, fitfun, sensitivities: list, fitness_dict: dict = {}, pop_id: str = '', print_progress: bool = True) -> (list, dict)
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.ga_param.Genetic_Algorithm.main

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.ga_param.Genetic_Algorithm.main

   .. py:method:: _check_if_kcats_are_in_bounds(individual) -> Individual
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.ga_param.Genetic_Algorithm._check_if_kcats_are_in_bounds

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.ga_param.Genetic_Algorithm._check_if_kcats_are_in_bounds

   .. py:method:: _get_best_individual_from_population(population: list)
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.ga_param.Genetic_Algorithm._get_best_individual_from_population

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.ga_param.Genetic_Algorithm._get_best_individual_from_population

   .. py:method:: _clone_elite(pop_to_clone: list, toolbox: deap.base.Toolbox) -> list
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.ga_param.Genetic_Algorithm._clone_elite

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.ga_param.Genetic_Algorithm._clone_elite
