#!/usr/bin/env python
# coding: utf-8
#usage: ./clonality_corr.py data_name data_name_extra data_root_dir var_name index result_name
#data_name: gse236581
#data_name_extra = 
#data_dir:/home/dreadcupper/jaewon/perturb_clonal_method/anaylsis/gse236581/
#var_name: var2000
#index: yes (include singlets) or no (only multiplets)
#result_name: cinema_result

import sys, os
import scanpy as sc
import pandas as pd
import numpy as np

import scipy.sparse as sp
from scipy.stats import spearmanr


data_name = sys.argv[1]
extra = sys.argv[2]
data_dir = sys.argv[3]
var_name = sys.argv[4]
index = sys.argv[5]
result_name = sys.argv[6]


tmp_result = pd.read_csv(data_dir + var_name + "/" + data_name + "_" + result_name, sep = "\t", index_col = 0)
pat_list = list(tmp_result.columns)
total_result = pd.DataFrame()
os.chdir(data_dir + var_name)


for pat in pat_list:
    print(pat)
    adata = sc.read_h5ad(data_dir + var_name + "/" + data_name + "_" + pat + "_" +  extra + '.h5ad')

    if index == "yes":
        adata1 = adata[(adata.obs["Tissue"] == "Tumor")]
    else:
        adata1 = adata[(adata.obs["Tissue"] == "Tumor") & (adata.obs["expand"] == "yes")]

    if adata1.n_obs < 50:
        total_result[pat] = pd.DataFrame(["NA","NA"])
        continue

    sc.pp.neighbors(adata1, use_rep = 'salient_rep', n_neighbors=5, knn = True)
    knn_graph = adata1.obsp['connectivities']
    clonality = adata1.obs['clonality'].values
    neighbor_sums = knn_graph@clonality
    neighbor_counts = np.array(knn_graph.sum(1)).flatten()
    mean_clonality = neighbor_sums / neighbor_counts

    corr, p_val = spearmanr(clonality, mean_clonality)
    total_result[pat] = pd.DataFrame([corr, p_val])

    ## random


total_result.index = ["clonality_corr", "clonality_pval"]
new_result = pd.concat([tmp_result, total_result])
new_result.to_csv(data_dir + var_name + "/" + data_name + "_" + result_name +"_corr", sep = "\t")

