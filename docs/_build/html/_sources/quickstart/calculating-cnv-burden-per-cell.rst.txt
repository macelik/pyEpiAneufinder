Calculating CNV Burden Per Cell
===============================

CNV burden is the aneuploidy score computed at the individual-cell level. It
can be used as an additional feature when distinguishing tumor cells.

.. code-block:: python

   import pandas as pd
   import pyEpiAneufinder as pea

   res = pd.read_csv(
       "results_sample_data/outs/result_table.tsv.gz",
       sep="\t",
       index_col=0,
   )
   cnv_burden = pea.compute_cnv_burden_cell(res)
