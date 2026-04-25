Quickstart
==========

Run the main workflow with a fragment file:

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
