import pandas as pd
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import seaborn as sns

import anndata as ad

def split_subclones(res, split_val, criterion="maxclust",
                    dist_metric="euclidean", linkage_method="ward"):
    """
    Split cells into subclones based on hierarchical clustering of CNV profiles.

    Parameters
    ----------
    res : pandas.DataFrame
        Result table from :func:`pyEpiAneufinder.epiAneufinder` containing
        ``seq``, ``start``, and ``end`` columns plus one CNV profile column per cell.
    split_val : int | float
        Threshold used to cut the clustering tree. With
        ``criterion="maxclust"``, this is the requested number of subclones.
        With ``criterion="distance"``, it is the clustering distance cutoff.
    criterion : {"maxclust", "distance"}, optional
        Criterion passed to :func:`scipy.cluster.hierarchy.fcluster`.
    dist_metric : str, optional
        Distance metric used when comparing cell-level CNV profiles.
    linkage_method : str, optional
        Linkage method used to build the hierarchical clustering tree.

    Returns
    -------
    pandas.DataFrame
        Data frame with two columns: ``barcode`` and ``subclone``.
    """

    # Remove position information (only CNVs kept)
    data_matrix = res.drop(columns=["seq","start","end"])

    # Calculate pairwise distances between cells
    dist_matrix = pdist(data_matrix.T, metric=dist_metric)

    # Normalize the distance by the number of bins (=> mean absolute error)
    fract_deviation = dist_matrix / res.shape[0]

    # Hierarchical clustering
    hc_cluster = linkage(fract_deviation, method=linkage_method)

    # Split into groups
    cl_members = fcluster(Z=hc_cluster, t=split_val, criterion=criterion)
    
    clones = pd.DataFrame({"barcode": data_matrix.columns.values,
                           "subclone": cl_members})

    return clones
