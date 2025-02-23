import math
import pandas as pd
import scanpy as sc
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
import time
import tracemalloc
from tensorflow.keras.models import Model
from tensorflow.keras.layers import (
    Input,
    Dense,
    BatchNormalization,
    Dropout,
    GaussianNoise,
)
from tensorflow.keras.regularizers import l2, l1
from tensorflow.keras.callbacks import LearningRateScheduler, EarlyStopping
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.losses import MeanSquaredLogarithmicError
from processing.pp import data_init, data_processing
from denoising.DAE import autoencoder
from clustering.lv import cluster_with_louvain
from marker.mk import CellTICS_marker_genes, get_unique_marker_genes
from analysis.ds import downstream_analysis


dataset_path = "../data/Usoskin.csv"

adata = data_init(dataset_path)
adata = data_processing(adata)


high_markers, low_markers = CellTICS_marker_genes(adata, thr1=0.0, thr2= 0.95)

high_markers, low_markers = get_unique_marker_genes(high_markers, low_markers)

adata = adata[:, adata.var_names.isin(high_markers)].copy()

expression_matrix = pd.read_csv("./pollen/count_table.csv", index_col=0).astype('float64')
expression_matrix.columns = [cell.split('.')[0] for cell in expression_matrix.columns]
expression_matrix = expression_matrix[expression_matrix.sum(axis = 1) > 0]

metadata = pd.read_csv("./pollen/metadata.csv", index_col=0)
metadata.set_index('Sample', inplace=True)
expression_matrix = expression_matrix.T

result = downstream_analysis(expression_matrix, metadata)