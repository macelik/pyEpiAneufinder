Input and Output
================

Inputs
------

- fragment file sorted by barcode and genomic position, or cellranger matrix directory
- reference genome FASTA
- blacklist BED file

Outputs
-------

- ``result_table.csv`` with CNV state per cell and bin
- ``count_matrix.h5ad`` with intermediate matrix outputs
- ``Karyogram.png`` when plotting is enabled
