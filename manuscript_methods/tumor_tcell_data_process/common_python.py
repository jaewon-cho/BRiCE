#!/usr/bin/env python
# coding: utf-8
#usage: ./common_python.py data_name data_name_extra data_root_dir var_name 
#data_name: gse236581
#data_name_extra = tcell_seurat
#data_dir:/home/dreadcupper/jaewon/perturb_clonal_method/anaylsis/gse236581/
#var_name: var2000


import sys, os
import scanpy as sc

data_name = sys.argv[1]
data_name_extra = sys.argv[2]
data_dir = sys.argv[3]
var_name = sys.argv[4]


adata = sc.read_h5ad(data_dir + data_name + "_" + data_name_extra +  '.h5ad')


if not os.path.exists(data_dir + "common_python/" + var_name):
    os.makedirs(data_dir + "common_python/" + var_name)
os.chdir(data_dir + "common_python/" + var_name)

pat_list = set(adata.obs["Patient"])

index = False
if var_name == "var_whole":
    index = True
else:
    var_gn = int(var_name[3:])

for pat in pat_list:
    print(pat)
    adata1 = adata[adata.obs["Patient"] == pat]
    sc.pp.normalize_total(adata1)
    sc.pp.log1p(adata1)
    if index:
        var_gn = adata.shape[1]
    sc.pp.highly_variable_genes(adata1, n_top_genes=var_gn)
    sc.tl.pca(adata1)
    sc.pp.neighbors(adata1, n_neighbors=10, n_pcs=30)
    sc.tl.umap(adata1)
    adata1.write("gse236581_" + pat + "_scapny.h5ad")


