Splitting Cells Into Clones
===========================

An approximate subclone assignment can be derived from the hierarchical
clustering using :func:`pyEpiAneufinder.split_subclones`.

.. code-block:: python

   import pandas as pd
   import pyEpiAneufinder as pea

   res = pd.read_csv(
       "results_sample_data/outs/result_table.tsv.gz",
       sep="\t",
       index_col=0,
   )
   clones = pea.split_subclones(res, split_val=4)

These clone labels can be shown as an annotation bar in
:func:`pyEpiAneufinder.karyo_gainloss`. The annotation table must be indexed by
cell barcode and must contain a column named ``annot``.

.. code-block:: python

   annot_dt = clones.copy()
   annot_dt.index = annot_dt.barcode
   annot_dt["annot"] = pd.Categorical("clone" + annot_dt.subclone.astype(str))

   pea.karyo_gainloss(
       res,
       outdir="results_sample_data/karyo_annot.png",
       title="Karyogram with annotation",
       annot_dt=annot_dt,
   )

.. image:: ../../sample_data/karyo_annot.png
   :alt: Karyogram with clone annotations
   :width: 900
