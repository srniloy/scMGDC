import math
import pandas as pd
import scanpy as sc
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from tensorflow.keras.regularizers import l2


# Define fphi - High markers
def fphi(x):

    # Step 1: Normalized by max value (x + 1) / (max(x) + 1)
    normalized = (x + 1).div(x.max(axis=1) + 1, axis=0)

    # Step 2: Suppression factor (zt)
    zt = (1 - normalized).sum(axis=1) / (x.shape[1] - 1)

    # Step 3: phi values after applying suppression factor
    zp = normalized.mul(zt, axis=0)

    return zp

# Define frho - Low markers
def frho(x):
    # Step 1: Normalized by min value (1 / (x + 1)) * (min(x) + 1)
    normalized = (1 / (x + 1)).mul(x.min(axis=1) + 1, axis=0)

    # Step 2: Enhancement factor (zt)
    zt = (1 - normalized).sum(axis=1) / (x.shape[1] - 1)

    # Step 3: rho values after applying enhancement factor
    zp = normalized.mul(zt, axis=0)

    return zp



# thr1 -> high marker quantile | thr2 -> Low marker quantile
def CellTICS_marker_genes(adata, thr1=0.95, thr2=0.9):
  print("\n====================== Identify Marker Gene With CellTICS ============================\n")
  print("Expression Values of Genes: \n")
  expr_data = pd.DataFrame(adata.X.T, index=adata.var_names, columns=adata.obs_names)
  louvain_labels = pd.Series(adata.obs['louvain'].values, index=adata.obs_names)
  print(expr_data.iloc[:10, :10])
  print(adata)

  print("\nAverage the cell expr data according to the clusters:\n")
  avg_expr_data = expr_data.T.groupby(louvain_labels).mean().T
  print(avg_expr_data.iloc[:10, :])

  phi = fphi(avg_expr_data)
  rho = frho(avg_expr_data)

  # Step 1: Extract Gene Names and Cluster Names --------------

  gnm = phi.index.values  # Gene names
  ctpnm = phi.columns.values  # Cluster names

  # Step 2: Convert DataFrames to Numpy Arrays ----------------

  phi = np.array(phi)
  rho = np.array(rho)

  # print("Final phi :")
  # print(phi[:10, :])
  # print("Final rho :")
  # print(rho[:10, :])

  # Step 3: Calculate Number of Marker Genes ------------------

  nummkg1 = math.ceil((1 - thr1) * phi.shape[0])  # High marker threshold
  nummkg2 = math.ceil((1 - thr2) * phi.shape[0])  # Low marker threshold


  # Step 4: Compute Quantiles for Each Cluster ----------------

  alpha = []  # Thresholds for phi
  beta = []  # Thresholds for rho
  for i in range(0, phi.shape[1]):
      alpha.append(np.quantile(phi[:, i], thr1))  # 95th percentile for phi
      beta.append(np.quantile(rho[:, i], thr2))  # 90th percentile for rho


  # Step 5: Identify High and Low Marker Genes ----------------

  mkh = []  # High markers
  mkl = []  # Low markers
  for i in range(0, phi.shape[1]):
      # print(f"alpha[{i}]:", alpha[i])
      high_markers = gnm[phi[:, i] >= alpha[i]][0:nummkg1]
      low_markers = gnm[rho[:, i] >= beta[i]][0:nummkg2]
      mkh = np.concatenate([mkh, high_markers], axis=0) if len(mkh) > 0 else high_markers
      mkl = np.concatenate([mkl, low_markers], axis=0) if len(mkl) > 0 else low_markers

  # Step 6: Reshape and Convert Results to DataFrames -------------

  mkh = mkh.reshape(nummkg1, phi.shape[1])
  mkl = mkl.reshape(nummkg2, rho.shape[1])
  mkh = pd.DataFrame(mkh, columns=ctpnm)
  mkl = pd.DataFrame(mkl, columns=ctpnm)

  print("\nHigh Marker Genes for each cluster/type: \n")
  print(mkh.iloc[:10, :])

  print("\nLow Marker Genes for each cluster/type: \n")
  print(mkl.iloc[:10, :])

  return mkh, mkl



def get_unique_marker_genes(high_markers, low_markers):
  unique_high_markers = high_markers
  unique_low_markers = low_markers
  unique_high_markers = pd.unique(unique_high_markers.values.ravel()).tolist()
  unique_low_markers = pd.unique(unique_low_markers.values.ravel()).tolist()
  print(f"\nTotal Unique High Marker Gene for all together {len(unique_high_markers)}: ")
  # print(unique_high_markers)
  print(f"\nTotal Unique Low Marker Gene for all together {len(unique_low_markers)}: ")
  # print(unique_low_markers)
  return unique_high_markers, unique_low_markers
