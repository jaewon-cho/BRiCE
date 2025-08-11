#!/usr/bin/env python
# coding: utf-8
#usage: ./performance_original.py folder_name index 

#folder_name 
#index: yes (include singlets) or no (only multiplets)

import sys, os
import pandas as pd
sys.path.append("/home/dreadcupper/jaewon/perturb_clonal_method/program/performance.py")
from performance import entropy_process, clonality_correlation
import glob
import scanpy as sc

index = sys.argv[2]
folder_name = sys.argv[1]

latent_category = "X_embedding"
category = "expand" 
k = 5


file_lst = sorted(glob.glob(os.path.join(folder_name, "*.h5ad")))
total_result =pd.DataFrame()

for file_name in file_lst:
    adata = sc.read_h5ad(file_name)
    pat = file_name.replace(".h5ad", "").split("_")[-2]
    print(pat)
    if index == "yes":
        obj = adata[(adata.obs["condition"] == "perturb")]
        obj1 = obj
    else:
        obj = adata[(adata.obs["condition"] == "perturb") & (adata.obs["expand"] == "yes")]
        obj1 = adata[(adata.obs["condition"] == "perturb")]
    
    
    ent_res = entropy_process(obj1, category)
    corr_res = clonality_correlation(obj, k)

    total_result[pat] = pd.DataFrame([ent_res.mean(), corr_res[0], corr_res[1]])

total_result.index = ["expand_entropy", "expand_corr", "corr_pval"]
total_result.to_csv(folder_name + "performance", sep = "\t")


