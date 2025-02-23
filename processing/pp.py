import scanpy as sc
import numpy as np


def data_init(path, type='Pollen'):
  print("\n==========================  Data Initialization  =================================\n")
  if type == 'Pollen' or 'Kolod':
    data = pd.read_csv(path, index_col=0).astype('float64')
    clusters = data.columns
    clusters = [int(x) for x in [float(x) for x in clusters]]
    adata = sc.AnnData(X=(data.values.T))
    adata.var_names = data.index
    adata.obs_names = [f"cell-{i}" for i in range(1, len(clusters) + 1)]
    adata.obs["cluster"] = clusters
    print("Loaded Dataset Structure: ")
    print((pd.DataFrame(adata.X, index=adata.obs.index, columns=adata.var.index)).iloc[:10, :10])
    print(adata)
    print(f'Total Clusters: \n{np.unique(adata.obs["cluster"])}')
    print(f'Clusters: \n{adata.obs["cluster"]}')
    sc.pp.pca(adata, n_comps=50) # Principal Component Analysis
    sc.pp.neighbors(adata, n_pcs=50) # Calculating Neighbors
    sc.tl.umap(adata) # Visualize the clusters
    print("Existing Cluster")
    sc.pl.umap(adata, color='cluster')
  return adata

def data_processing(adata):
  print("\n==========================  Data Processing  =================================\n")
  sc.pp.filter_cells(adata, min_genes=200)
  sc.pp.filter_genes(adata, min_cells=3)
  sc.pp.normalize_per_cell(adata)
  sc.pp.log1p(adata, copy=True)
  return adata