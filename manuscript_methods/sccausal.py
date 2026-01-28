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
from scCausalVI import scCausalVIModel
import torch

import seaborn as sns
import matplotlib.pyplot as plt

from math import log
import random
import bisect
from scipy.stats import wilcoxon
from scipy.stats import spearmanr

scvi.settings.seed = 0


sample_list_file = sys.argv[1]
data_name = sys.argv[2]
data_root_dir = sys.argv[3]
var_name = sys.argv[4]
save_root_dir = sys.argv[5]
perturb_only = sys.argv[6]
expand_only = sys.argv[7]


# In[2]:


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


# In[3]:


def entropy_process(latent_obj, category):
    #category: expand or Tissue
    dataset=latent_obj
    sc.pp.neighbors(dataset, use_rep='latent_t', n_neighbors=31, knn = True)
    #n_neighbor: batchbench
    
    knn_graph = dataset.obsp['distances']
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


# In[4]:


def entropy_rand_process(dataset, category):
    #dataset=latent_obj
    #
    #n_neighbor: batchbench
    
    knn_graph = dataset.obsp['distances']
    #transforming csr_matrix to dataframe
    df = pd.DataFrame(knn_graph.toarray())
    rand_ind =list(range(len(df)))
    random.shuffle(rand_ind)
    
    kwargs = {}
    #batch vector(batch id of each cell)
    kwargs['batch_vector'] = pd.Series(dataset.obs[category], index = dataset.obs[category].index[rand_ind])
    #modify index of batch vector so it coincides with matrix's index
    kwargs['batch_vector'].index = range(0,len(kwargs['batch_vector']))
    #number of batches
    kwargs['N_batches'] = len(dataset.obs[category].astype('category').cat.categories)
    entropy_result = df.apply(shannon_entropy, axis=0, args=(kwargs['batch_vector'],kwargs['N_batches']))
    
    return(entropy_result)

def clonality_correlation(adata1, k, latent = None):
    sc.pp.neighbors(adata1, use_rep = 'latent_t', n_neighbors=k+1, knn = True)
    knn_binary = (adata1.obsp["distances"] > 0).astype(int)
    clonality = adata1.obs['clonality'].values
    neighbor_sums = knn_binary@clonality
    neighbor_counts = np.array(knn_binary.sum(1)).flatten()
    mean_clonality = neighbor_sums / neighbor_counts

    corr, p_val = spearmanr(clonality, mean_clonality)
    return([corr, p_val])


def sccausal_process(pat):
    if os.path.isfile(data_root_dir + var_name + "/" + data_name + '_' + pat + '_scapny.h5ad'):
        adata = sc.read_h5ad(data_root_dir + var_name + "/" + data_name + '_' + pat + '_scapny.h5ad')
    else:
        return "apple"

    k = 5

    if perturb_only == "yes":
        condition_key = 'expand'
        conditions = ['no', 'yes',]
        control_key='no'
        group_indices_list = [np.where(adata.obs[condition_key]==group)[0] for group in conditions]
    else:
        condition_key = 'condition'
        conditions = ['unperturb', 'perturb',]
        control_key='unperturb'
        group_indices_list = [np.where(adata.obs[condition_key]==group)[0] for group in conditions]

    #adata.layers["counts"] = adata.X.copy()
    scCausalVIModel.setup_anndata(adata, condition_key=condition_key, layer='counts',)

# For multi-batch data, register batch index of batch_key in obs:
# scCausalVIModel.setup_anndata(adata, batch_key='stim', condition_key=condition_key, layer='counts',)

    condition2int = adata.obs.groupby(condition_key, observed=False)['_scvi_condition'].first().to_dict()

    model = scCausalVIModel(
    adata,
    condition2int=condition2int,
    control=control_key,
    n_background_latent=10,
    n_te_latent=10,
    n_layers=2,
    n_hidden=128,
    use_mmd=True,
    mmd_weight=10,
    norm_weight=0.2,
    )


    model.train(
    group_indices_list,
    use_gpu=use_gpu,
    max_epochs=200,
    )
   
    adata.obsm['latent_bg'], adata.obsm['latent_t'] = model.get_latent_representation()

    adata.write(data_name + "_" + pat + "_sccausal.h5ad")


    if perturb_only == "yes":
        latent_after_entropy = entropy_process(adata,'expand')

        new_obj1 = adata.copy()
        expand_entropy = entropy_process(new_obj1, "expand")
        if expand_only == "yes":
            new_obj11 = new_obj1[new_obj1.obs["expand"] == "yes"]
        else:
            new_obj11 = new_obj1.copy()


    else:
        latent_after_entropy = entropy_process(adata,'condition')

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

    results = [latent_after_entropy.mean(), expand_entropy.mean(), corr_res[0], corr_res[1]]
    return(results)


###
use_gpu = torch.cuda.is_available()





###
pat_list = list(line.strip() for line in open(sample_list_file))

if not os.path.exists(save_root_dir + var_name):
    os.makedirs(save_root_dir + var_name)
os.chdir(save_root_dir + var_name)
total_result = pd.DataFrame()

for pat in pat_list:
    print(pat)
    pat_res = sccausal_process(pat)
    if pat_res == "apple":
        continue
    #sinkhorn_knopp does not fit

    df = pd.DataFrame(pat_res)
    df.index = ["latent_after_entropy","expand_entropy","expand_corr", "corr_pval"]
    total_result[pat] = pd.Series(pat_res)
total_result.index = ["latent_after_entropy","expand_entropy","expand_corr", "corr_pval"]
total_result.to_csv(data_name + "_sccausal_result", sep = "\t")



