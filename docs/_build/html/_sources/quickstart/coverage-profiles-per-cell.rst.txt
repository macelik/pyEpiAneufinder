Coverage Profiles Per Cell
==========================

Individual-cell CNV calls can be explored further with
:func:`pyEpiAneufinder.plot_single_cell_profile`.

.. code-block:: python

   import pyEpiAneufinder as pea

   pea.plot_single_cell_profile(
       outdir="results_sample_data",
       cell_name="AGTCCGGTCCACACCT-1",
       plot_path="results_sample_data/somy_profile_cell.png",
   )

.. image:: ../../sample_data/somy_profile_cell.png
   :alt: Per-cell coverage profile from the example dataset
   :width: 900
