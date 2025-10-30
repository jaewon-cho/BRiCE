#!/usr/bin/env python
# coding: utf-8

# In[ ]:

#usage: sccausal_run.py sample_list_file data_name data_root_dir var_name save_root_dir perturb_only expand_only
#sample_file_list:
#data_name: gse236581
# data_root_dir include last  /
# /home/dreadcupper/jaewon/perturb_clonal_method/anaylsis/gse236581/common_python/
#var_name: var2000
# save_root_dir include last  /
# /home/dreadcupper/jaewon/perturb_clonal_method/anaylsis/gse236581/sccausal/

#perturb_only: yes or no
#yes: perturb only
#expand_only: yes or no
#yes: expand only for measuring clonality correlation

import os,sys
import scanpy as sc
import scvi
import anndata as ad
import numpy as np
import pandas as pd
import torch

import matplotlib.pyplot as plt

from math import log
import random
import bisect
from scipy.stats import wilcoxon
from scipy.stats import spearmanr

from scDisInFact import scdisinfact, create_scdisinfact_dataset
from scDisInFact import utils

#from umap import UMAP
from sklearn.decomposition import PCA
import scipy.sparse as sp

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

sample_list_file = sys.argv[1]
data_name = sys.argv[2]
data_root_dir = sys.argv[3]
var_name = sys.argv[4]
save_root_dir = sys.argv[5]
perturb_only = sys.argv[6]
expand_only = sys.argv[7]

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
    sc.pp.neighbors(dataset,use_rep= "latent", n_neighbors=30, knn = True)
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
    sc.pp.neighbors(adata1, use_rep= "latent",n_neighbors=k, knn = True)
    knn_graph = adata1.obsp['connectivities']
    clonality = adata1.obs['clonality'].values
    neighbor_sums = knn_graph@clonality
    neighbor_counts = np.array(knn_graph.sum(1)).flatten()
    mean_clonality = neighbor_sums / neighbor_counts

    corr, p_val = spearmanr(clonality, mean_clonality)
    return([corr, p_val])

def scdisinfact_process(pat):

    reg_mmd_comm = 1e-4
    reg_mmd_diff = 1e-4
    reg_kl_comm = 1e-5
    reg_kl_diff = 1e-2
    reg_class = 1
    reg_gl = 1

    #Ks = [8, 2, 2]
    Ks = [1,2]

    batch_size = 64
    nepochs = 100
    interval = 10
    lr = 5e-4
    lambs = [reg_mmd_comm, reg_mmd_diff, reg_kl_comm, reg_kl_diff, reg_class, reg_gl]

    if os.path.isfile(data_root_dir + var_name + "/" + data_name + '_' + pat + '_scapny.h5ad'):
        adata = sc.read_h5ad(data_root_dir + var_name + "/" + data_name + '_' + pat + '_scapny.h5ad')
    else:
        return "apple"

    k = 5

    counts = adata.layers['counts'].copy()
    meta_cells = adata.obs.copy()


    if perturb_only == "yes":
        data_dict = create_scdisinfact_dataset(counts, meta_cells, condition_key = ["expand"], batch_key = "individual")
    else:
        data_dict = create_scdisinfact_dataset(counts, meta_cells, condition_key = ["condition"], batch_key = "individual")
    
    model = scdisinfact(data_dict = data_dict, Ks = Ks, batch_size = batch_size, interval = interval, lr = lr, reg_mmd_comm = reg_mmd_comm, reg_mmd_diff = reg_mmd_diff, reg_gl = reg_gl, reg_class = reg_class, reg_kl_comm = reg_kl_comm, reg_kl_diff = reg_kl_diff, seed = 0, device = device)
    model.train()
    losses = model.train_model(nepochs = nepochs, recon_loss = "NB")

    pert = data_dict["datasets"][0]
    unpert = data_dict["datasets"][1]


    dict_inf1 = model.inference(counts = pert.counts_norm.to(model.device), batch_ids = pert.batch_id[:,None].to(model.device), print_stat = True)
    dict_inf2 = model.inference(counts = unpert.counts_norm.to(model.device), batch_ids = unpert.batch_id[:,None].to(model.device), print_stat = True)

    z_d1 = dict_inf1["mu_d"]
    z_d2 = dict_inf2["mu_d"]
    pca_op = PCA(n_components = 2)

    z_np1 = z_d1[0].detach().cpu().numpy()
    z_np2 = z_d2[0].detach().cpu().numpy()
    latent = pca_op.fit_transform(np.concatenate([z_np1, z_np2], axis = 0))
    #latent = pca_op.fit_transform(z_np)
    cells1 = data_dict["meta_cells"][0].index
    cells2 = data_dict["meta_cells"][1].index

    all_cells = pd.Index(cells1.tolist() + cells2.tolist(), name="cell_id")

    latent_df = pd.DataFrame(latent, index=all_cells)
    latent_df = latent_df.loc[adata.obs_names]
    adata.obsm["latent"] = latent_df.to_numpy()

    adata.write(data_name + "_" + pat + "_scdisinfact.h5ad")

    if perturb_only == "yes":
        new_obj1 = adata.copy()
        expand_entropy = entropy_process(new_obj1, "expand")
        if expand_only == "yes":
            new_obj11 = new_obj1[new_obj1.obs["expand"] == "yes"]
        else:
            new_obj11 = new_obj1.copy()


    else:
        new_obj1 = adata[adata.obs["condition"] == "perturb"]
        expand_entropy = entropy_process(new_obj1, "expand")
        if expand_only == "yes":
            new_obj11 = new_obj1[(new_obj1.obs["expand"].values == "yes") & (new_obj1.obs["condition"].values == "perturb")]
        else:
            new_obj11 = new_obj1[new_obj1.obs["condition"] == "perturb"]

    if new_obj11.n_obs < 50:
        corr_res = ["NA", "NA"]
    else:
        corr_res = clonality_correlation(new_obj11, k)

    results = [expand_entropy.mean(), corr_res[0], corr_res[1]]
    return(results)





pat_list = list(line.strip() for line in open(sample_list_file))

if not os.path.exists(save_root_dir + var_name):
    os.makedirs(save_root_dir + var_name)
os.chdir(save_root_dir + var_name)
total_result = pd.DataFrame()

for pat in pat_list:
    print(pat)
    pat_res = scdisinfact_process(pat)
    if pat_res == "apple":
        continue
    #sinkhorn_knopp does not fit

    df = pd.DataFrame(pat_res)
    df.index = ["expand_entropy","expand_corr", "corr_pval"]
    total_result[pat] = pd.Series(pat_res)
total_result.index = ["expand_entropy","expand_corr", "corr_pval"]
total_result.to_csv(data_name + "_scdisinfact_result", sep = "\t")



