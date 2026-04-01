:py:mod:`PAMparametrizer.utils.pamparametrizer_visualization`
=============================================================

.. py:module:: PAMparametrizer.utils.pamparametrizer_visualization

.. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_visualization
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`plot_simulation <PAMparametrizer.utils.pamparametrizer_visualization.plot_simulation>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_visualization.plot_simulation
          :summary:
   * - :py:obj:`plot_flux_vs_experiment <PAMparametrizer.utils.pamparametrizer_visualization.plot_flux_vs_experiment>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_visualization.plot_flux_vs_experiment
          :summary:
   * - :py:obj:`plot_valid_data <PAMparametrizer.utils.pamparametrizer_visualization.plot_valid_data>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_visualization.plot_valid_data
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`FONTSIZE <PAMparametrizer.utils.pamparametrizer_visualization.FONTSIZE>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_visualization.FONTSIZE
          :summary:
   * - :py:obj:`RXN_NAME_MAPPER <PAMparametrizer.utils.pamparametrizer_visualization.RXN_NAME_MAPPER>`
     - .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_visualization.RXN_NAME_MAPPER
          :summary:

API
~~~

.. py:data:: FONTSIZE
   :canonical: PAMparametrizer.utils.pamparametrizer_visualization.FONTSIZE
   :value: 16

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_visualization.FONTSIZE

.. py:data:: RXN_NAME_MAPPER
   :canonical: PAMparametrizer.utils.pamparametrizer_visualization.RXN_NAME_MAPPER
   :value: None

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_visualization.RXN_NAME_MAPPER

.. py:function:: plot_simulation(fig, axs, fluxes: pandas.DataFrame, substrate_rates: list, reactions_to_plot: list, iteration: int = 0, color: int = None, max_iteration: int = 2, label: str = None, return_color: bool = False, plotting_kwargs: typing.Dict[str, typing.Any] = {}) -> matplotlib.pyplot.Figure
   :canonical: PAMparametrizer.utils.pamparametrizer_visualization.plot_simulation

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_visualization.plot_simulation

.. py:function:: plot_flux_vs_experiment(ax, parametrizer, color, fontsize)
   :canonical: PAMparametrizer.utils.pamparametrizer_visualization.plot_flux_vs_experiment

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_visualization.plot_flux_vs_experiment

.. py:function:: plot_valid_data(parametrizer, axs=None, fig=None, fontsize: int = 12)
   :canonical: PAMparametrizer.utils.pamparametrizer_visualization.plot_valid_data

   .. autodoc2-docstring:: PAMparametrizer.utils.pamparametrizer_visualization.plot_valid_data
