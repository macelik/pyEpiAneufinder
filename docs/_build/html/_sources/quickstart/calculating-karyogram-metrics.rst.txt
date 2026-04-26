Calculating Karyogram Metrics
=============================

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

.. image:: ../../sample_data/scatter_aneu_heterogen.png
   :alt: Scatter plot of chromosome-level aneuploidy versus heterogeneity
   :width: 700
