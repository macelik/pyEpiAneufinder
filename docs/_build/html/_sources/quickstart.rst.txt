Quickstart
==========

This page walks through the main ``pyEpiAneufinder`` workflow using the bundled
example dataset and the primary public helper functions documented for the first
release.

Executing the Program
---------------------

The main workflow is executed through :func:`pyEpiAneufinder.epiAneufinder`.
When using a fragment file as input, the file should be ordered first by cell
barcode and then by genomic position. ``pyEpiAneufinder`` can do this sorting
internally when ``sort_fragment=True``.

Alternatively, the ``fragment_file`` argument can point to a Cell Ranger-style
matrix directory containing the peaks, barcodes, and matrix files. In that
case, set ``cellRangerInput=True``. Fragment-file input is still the preferred
path when available.

By default, the workflow performs GC correction. The correction can be disabled
with ``GC=False``, but the default is recommended.

The workflow also needs:

- a reference genome FASTA file
- a blacklist BED file
- a window size
- an output directory
- optional chromosomes to exclude
- a minimum fragment cutoff for filtering cells

For example:

.. code-block:: python

   import pyEpiAneufinder as pea

   pea.epiAneufinder(
       fragment_file="sample_data/sample.tsv.gz",
       outdir="results_sample_data",
       genome_file="hg38.fa.gz",
       blacklist="sample_data/hg38-blacklist.v2.bed",
       windowSize=100000,
       ncores=1,
       exclude=["chrX", "chrY"],
       minFrags=20000,
       resume=False,
       cellRangerInput=False,
       GC=True,
       sort_fragment=True,
       remove_barcodes=None,
   )

The test dataset used in the examples is available in the repository under
``sample_data/``.

The main output is written to ``results_sample_data/outs/result_table.tsv.gz``.
This table contains the consensus CNV state per cell and genomic bin, encoded
as ``0=loss``, ``0.5=weak loss``, ``1=base``, ``1.5=weak gain``, and ``2=gain``.
Additional intermediate outputs are described on the :doc:`input-output` page.

The example run produces the following karyogram:

.. image:: ../sample_data/Karyogram.png
   :alt: Karyogram generated from the example dataset
   :width: 900

Splitting Cells Into Clones
---------------------------

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

.. image:: ../sample_data/karyo_annot.png
   :alt: Karyogram with clone annotations
   :width: 900

Calculating Karyogram Metrics
-----------------------------

The CNV profiles can be summarized with aneuploidy and heterogeneity metrics.
These are available for the full sample and per chromosome.

.. code-block:: python

   import matplotlib.pyplot as plt
   import pandas as pd
   import pyEpiAneufinder as pea
   import seaborn as sns

   res = pd.read_csv(
       "results_sample_data/outs/result_table.tsv.gz",
       sep="\t",
       index_col=0,
   )

   pea.compute_aneuploidy_across_sample(res)
   pea.compute_heterogeneity_across_sample(res)

   aneu_chrom = pea.compute_aneuploidy_by_chr(res)
   heterogen_chrom = pea.compute_heterogeneity_by_chr(res)

   plot_data = pd.DataFrame(
       {
           "chrom": aneu_chrom.columns.values,
           "aneu": aneu_chrom.iloc[0],
           "heterogen": heterogen_chrom.iloc[0],
       }
   )

   sns.scatterplot(x="aneu", y="heterogen", data=plot_data)

   for i in range(len(plot_data)):
       plt.annotate(
           plot_data["chrom"][i],
           (plot_data["aneu"][i], plot_data["heterogen"][i]),
       )

   plt.xlabel("Aneuploidy per chromosome")
   plt.ylabel("Heterogeneity per chromosome")
   plt.show()

For the example data, the resulting scatter plot looks like this:

.. image:: ../sample_data/scatter_aneu_heterogen.png
   :alt: Scatter plot of chromosome-level aneuploidy versus heterogeneity
   :width: 700

Calculating CNV Burden Per Cell
-------------------------------

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

Coverage Profiles Per Cell
--------------------------

Individual-cell CNV calls can be explored further with
:func:`pyEpiAneufinder.plot_single_cell_profile`.

.. code-block:: python

   import pyEpiAneufinder as pea

   pea.plot_single_cell_profile(
       outdir="results_sample_data",
       cell_name="AGTCCGGTCCACACCT-1",
       plot_path="results_sample_data/somy_profile_cell.png",
   )

.. image:: ../sample_data/somy_profile_cell.png
   :alt: Per-cell coverage profile from the example dataset
   :width: 900
