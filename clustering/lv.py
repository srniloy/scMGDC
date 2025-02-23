import numpy as np
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
import math
import pandas as pd
import scanpy as sc

def cluster_with_louvain(xdata, epochs, resolution=2.5):
  print("\n========================== Clustering With Louvain =================================\n")
  if "X_pca" not in xdata.obsm and "neighbors" not in xdata.uns:
    sc.pp.pca(xdata, n_comps=50)
    sc.pp.neighbors(xdata, n_pcs=50)

  for i in range(0, epochs):
    sc.tl.louvain(xdata, resolution)

    print(f"\nClustering epoch {i+1} ------------------")
    ari = adjusted_rand_score(xdata.obs['cluster'], xdata.obs['louvain'])
    print(f"Adjusted Rand Index (ARI): {ari:.4f}")

    # Compute Normalized Mutual Information (NMI)
    nmi = normalized_mutual_info_score(xdata.obs['cluster'], xdata.obs['louvain'])
    print(f"Normalized Mutual Information (NMI): {nmi:.4f}\n")
  sc.tl.umap(xdata)
  sc.pl.umap(xdata, color='louvain')
  return xdata