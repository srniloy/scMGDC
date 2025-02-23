import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


def downstream_analysis(expression_matrix, metadata):
    dds = DeseqDataSet(
        counts=expression_matrix, 
        metadata=metadata,  
        design_factors=["Condition"]
    )
    dds.deseq2()
    stat_res = DeseqStats(dds)
    stat_res.summary()
    stat_res
    results_df = stat_res.results_df
    sort_sg = significant_genes.dropna().sort_values('padj', ascending=False)
    pos_sg = sort_sg[sort_sg['log2FoldChange'] >= 0]
    neg_sg = sort_sg[sort_sg['log2FoldChange'] < 0]
    return (sort_sg, pos_sg, neg_sg)