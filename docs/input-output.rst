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

- ``outdir/outs/result_table.tsv.gz`` with the consensus Holmes+Watson CNV table
- ``outdir/count_matrix.h5ad`` with the filtered count matrix used for inference
- ``outdir/parameter_configuration.txt`` with the run configuration
- ``outdir/binned_genome.csv`` with the genome bins and GC annotations
- ``outdir/breakpoints.csv`` with the detected breakpoints per cell
- ``outdir/clusters.json`` with the inferred segment clusters per cell
- ``outdir/outs/Karyogram.png`` when plotting is enabled
- ``outdir/outs/integer_states.tsv.gz`` with integer CNV states
- ``outdir/outs/continuous_scores.tsv.gz`` with continuous CNV scores
- ``outdir/outs/result_table_holmes.tsv.gz`` with Holmes-mapped CNV states
- ``outdir/outs/result_table_watson.tsv.gz`` with Watson-mapped CNV states
- ``outdir/outs/scaling_factors.tsv.gz`` with per-cell scaling factors
