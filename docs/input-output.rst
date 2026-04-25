Input and Output
================

Inputs
------

- fragment file or Cell Ranger matrix directory; when using a fragment file, the
  workflow sorts it internally by barcode and genomic position by default
- reference genome FASTA
- blacklist BED file

Outputs
-------

- ``outdir/outs/result_table.tsv.gz`` with CNV state per cell and bin
- ``count_matrix.h5ad`` with intermediate matrix outputs
- ``outdir/outs/Karyogram.png`` when plotting is enabled
