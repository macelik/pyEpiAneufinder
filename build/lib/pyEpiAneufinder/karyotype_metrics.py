import pandas as pd
import numpy as np
from natsort import natsorted

def compute_cnv_burden_cell(df,offset=3):
    """
    Compute the fraction of altered bins for each cell.

    Parameters
    ----------
    df : pandas.DataFrame
        Result table from :func:`pyEpiAneufinder.epiAneufinder`.
    offset : int, optional
        Number of leading metadata columns to skip before the per-cell CNV matrix.

    Returns
    ------
    pandas.DataFrame
        Data frame with ``barcodes`` and ``cnv_burden`` columns, where
        ``cnv_burden`` is the fraction of bins whose state differs from 1.
    """
    cn_matrix = df.iloc[:, offset:].to_numpy()
    cnv_burden = pd.DataFrame({"barcodes": df.columns.values[offset:],
                "cnv_burden":np.mean(cn_matrix != 1, axis=0)})
    return cnv_burden

def compute_aneuploidy_across_sample(df, offset=3):
    """
    Compute the mean fraction of altered bins across the full sample.

    Parameters
    ----------
    df : pandas.DataFrame
        Result table from :func:`pyEpiAneufinder.epiAneufinder`.
    offset : int, optional
        Number of leading metadata columns to skip before the per-cell CNV matrix.

    Returns
    ------
    float
        Sample-level aneuploidy score computed as the fraction of entries not
        equal to 1.
    """
    cn_matrix = df.iloc[:, offset:].to_numpy()
    return np.mean(cn_matrix != 1)

def compute_aneuploidy_by_chr(df, offset=3):
    """
    Compute an aneuploidy score for each chromosome.

    Parameters
    ----------
    df : pandas.DataFrame
        Result table from :func:`pyEpiAneufinder.epiAneufinder`.
    offset : int, optional
        Number of leading metadata columns to skip before the per-cell CNV matrix.

    Returns
    ------
    pandas.DataFrame
        Single-row data frame with one column per chromosome.
    """
    cn_matrix = df.iloc[:, offset:].to_numpy()
    chroms = df['seq'].values
    unique_chroms = natsorted(np.unique(chroms))

    result = {}
    for chr_ in unique_chroms:
        idx = np.where(chroms == chr_)[0]
        result[chr_] = np.mean(cn_matrix[idx] != 1)

    df_result = pd.DataFrame(result, index=[0])  # index is just zero
    return df_result

#Help function for compute_heterogeneity_across_sample()
def compute_heterogeneity_array(arr):
    heterogeneity = np.zeros(arr.shape[0])
    for i, row in enumerate(arr):
        vals, counts = np.unique(row, return_counts=True)
        counts = np.sort(counts)[::-1]
        weights = np.arange(len(counts))
        heterogeneity[i] = np.sum(counts * weights) / arr.shape[1]
    return heterogeneity

def compute_heterogeneity_across_sample(df, offset=3):
    """
    Compute the mean per-bin heterogeneity score across the full sample.

    Parameters
    ----------
    df : pandas.DataFrame
        Result table from :func:`pyEpiAneufinder.epiAneufinder`.
    offset : int, optional
        Number of leading metadata columns to skip before the per-cell CNV matrix.

    Returns
    ------
    float
        Mean heterogeneity score across all genomic bins.
    """
    cn_matrix = df.iloc[:, offset:].to_numpy()
    heterogeneity = compute_heterogeneity_array(cn_matrix)
    return np.mean(heterogeneity)

def compute_heterogeneity_by_chr(df, offset=3):
    """
    Compute a heterogeneity score for each chromosome.

    Parameters
    ----------
    df : pandas.DataFrame
        Result table from :func:`pyEpiAneufinder.epiAneufinder`.
    offset : int, optional
        Number of leading metadata columns to skip before the per-cell CNV matrix.

    Returns
    ------
    pandas.DataFrame
        Single-row data frame with one column per chromosome.
    """
    cn_matrix = df.iloc[:, offset:].to_numpy()
    chroms = df['seq'].values
    unique_chroms = natsorted(np.unique(chroms))

    result = {}
    for chr_ in unique_chroms:
        idx = np.where(chroms == chr_)[0]
        heterogeneity = compute_heterogeneity_array(cn_matrix[idx])
        result[chr_] = np.mean(heterogeneity)

    df_result = pd.DataFrame(result, index=[0])  # index is just zero
    return df_result

