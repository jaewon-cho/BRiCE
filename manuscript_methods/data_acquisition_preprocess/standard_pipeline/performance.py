#!/usr/bin/env python
# coding: utf-8
#usage: ./performance.py 

import sys, os
import scanpy as sc
import pandas as pd
import numpy as np

import scipy.sparse as sp
from scipy.stats import spearmanr
import argparse
from math import log
import random
import bisect
from scipy.stats import wilcoxon

def shannon_entropy (x, b_vec, N_b):
    #x: corresponding cluster (if >0: it is a given cluster)
    #b_vec: Tissue
    #N_b: the number of Tissue category

    tabled_values = b_vec[x > 0].value_counts()/ len(b_vec[x >0]) #class 'pandas.core.series.Series'

    tabled_val = tabled_values.tolist()

    entropy = 0.0
    for element in tabled_val:
        if element != 0:
            entropy += element * log(element)

    entropy /= log(N_b)

    return(-entropy) #the entropy formula is the -sum, this is why we include the minus sign

def entropy_process(latent_obj, category):
    #category: expand or Tissue
    dataset=latent_obj
    sc.pp.neighbors(dataset, n_neighbors=30, knn = True)
    #n_neighbor: batchbench

    knn_graph = dataset.obsp['connectivities']
    #transforming csr_matrix to dataframe
    df = pd.DataFrame(knn_graph.toarray())
    kwargs = {}
    #batch vector(batch id of each cell)
    kwargs['batch_vector'] = dataset.obs[category]
    #modify index of batch vector so it coincides with matrix's index
    kwargs['batch_vector'].index = range(0,len(kwargs['batch_vector']))
    #number of batches
    kwargs['N_batches'] = len(dataset.obs[category].astype('category').cat.categories)
    entropy_result = df.apply(shannon_entropy, axis=0, args=(kwargs['batch_vector'],kwargs['N_batches']))

    return(entropy_result)



def clonality_correlation(adata1, k):
    sc.pp.neighbors(adata1, n_neighbors=k, knn = True)
    knn_graph = adata1.obsp['connectivities']
    clonality = adata1.obs['clonality'].values
    neighbor_sums = knn_graph@clonality
    neighbor_counts = np.array(knn_graph.sum(1)).flatten()
    mean_clonality = neighbor_sums / neighbor_counts

    corr, p_val = spearmanr(clonality, mean_clonality)
    return([corr, p_val])



