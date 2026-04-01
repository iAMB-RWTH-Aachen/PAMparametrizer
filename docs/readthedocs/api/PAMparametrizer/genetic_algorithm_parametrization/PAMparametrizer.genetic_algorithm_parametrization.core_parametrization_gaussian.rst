:py:mod:`PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian`
=========================================================================================

.. py:module:: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian

.. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian
   :allowtitles:

Module Contents
---------------

Classes
~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`GAPO <PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO>`
     - .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO
          :summary:
   * - :py:obj:`MyFitness <PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.MyFitness>`
     -

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`print_time <PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.print_time>`
     - .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.print_time
          :summary:

API
~~~

.. py:data:: print_time
   :canonical: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.print_time
   :value: None

   .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.print_time

.. py:class:: GAPO(model=None, enzymes_to_eval: dict = {}, valid_data: typing.Union[cobra.DictList, dict] = {}, kcat_constraints_table: pandas.DataFrame = KcatConstraintConfigTable(), sector_configs_per_substrate: dict = {}, fitness_class='Fitfun_params_uniform', mutation_probability=0.5, mutation_rate=0.05, population_size=30, crossover_probability=0.8, number_generations=20, number_gene_flow_events=10, processes=2, time_limit=600, init_attribute_probability=0, fixed_attributes=[], folderpath_save=Path('Results'), filename_save='ga_results', overwrite_intermediate_results=True, objective_id='BIOMASS', sigma_denominator: int = 10, substrate_uptake_rates={'EX_glc__D_e': [0.7, 11.3]}, substrate_uptake_id='EX_glc__D_e', error_weights: dict = {}, print_progress=True)
   :canonical: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO

   .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO

   .. rubric:: Initialization

   .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO.__init__

   .. py:method:: start() -> None
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO.start

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO.start

   .. py:method:: restart(filepath_previous_pop: typing.Union[list, str])
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO.restart

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO.restart

   .. py:method:: restart_with_different_individuals(previous_enzyme_values: list[dict]) -> None
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO.restart_with_different_individuals

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO.restart_with_different_individuals

   .. py:method:: _init_deap_toolbox()
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO._init_deap_toolbox

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO._init_deap_toolbox

   .. py:method:: _init_deap_toolbox_mutation(toolbox)
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO._init_deap_toolbox_mutation

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO._init_deap_toolbox_mutation

   .. py:method:: _init_deap_individual()
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO._init_deap_individual

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO._init_deap_individual

   .. py:method:: _init_deap_fitness()
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO._init_deap_fitness

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO._init_deap_fitness

   .. py:method:: _copy_deap_individual(toolbox, to_copy)
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO._copy_deap_individual

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO._copy_deap_individual

   .. py:method:: _parallel_gene_flow(pops, toolbox, start_time: float, previous_drifts: float = 0) -> typing.List[list]
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO._parallel_gene_flow

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO._parallel_gene_flow

   .. py:method:: _save_population(pop, suffix: str = '') -> None
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO._save_population

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.GAPO._save_population

.. py:class:: MyFitness(values=())
   :canonical: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.MyFitness

   Bases: :py:obj:`deap.base.Fitness`

   .. py:attribute:: weights
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.MyFitness.weights
      :value: (1,)

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.MyFitness.weights

   .. py:method:: _wsum()
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.MyFitness._wsum

      .. autodoc2-docstring:: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.MyFitness._wsum

   .. py:method:: __le__(other)
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.MyFitness.__le__

   .. py:method:: __lt__(other)
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.MyFitness.__lt__

   .. py:method:: __eq__(other)
      :canonical: PAMparametrizer.genetic_algorithm_parametrization.core_parametrization_gaussian.MyFitness.__eq__
