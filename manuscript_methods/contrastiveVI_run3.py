#!/usr/bin/env python
# coding: utf-8

# In[1]:
#usage: ./contrastiveVI_run.py sample_list_file data_name data_root_dir var_name save_root_dir perturb_only expand_only
#sample_file_list:
#data_name: gse236581
# data_root_dir include last  /
# /home/dreadcupper/jaewon/perturb_clonal_method/anaylsis/gse236581/common_python/
#var_name: var2000
# save_root_dir include last  /
# /home/dreadcupper/jaewon/perturb_clonal_method/analysis/contrastiveVI/

#perturb_only: yes or no
#yes: perturb only
#expand_only: yes or no
#yes: expand only for measuring clonality correlation


from contrastive_vi.model import ContrastiveVI
from pytorch_lightning.utilities.seed import seed_everything
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

import argparse
from math import log
import random
import bisect
from scipy.stats import wilcoxon
import pandas as pd
import sys,os
from scipy.stats import spearmanr


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


def entropy_process(latent_obj, latent_category, category):
    #category: expand or Tissue
    dataset=latent_obj
    sc.pp.neighbors(dataset, use_rep= latent_category, n_neighbors=31, knn = True)
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
    sc.pp.neighbors(adata1, use_rep = 'salient_rep', n_neighbors=k+1, knn = True)
    knn_binary = (adata1.obsp["distances"] > 0).astype(int)
    clonality = adata1.obs['clonality'].values
    neighbor_sums = knn_binary@clonality
    neighbor_counts = np.array(knn_binary.sum(1)).flatten()
    mean_clonality = neighbor_sums / neighbor_counts

    corr, p_val = spearmanr(clonality, mean_clonality)
    return([corr, p_val])


# In[5]:


def contrastive_process(pat):
    if os.path.isfile(data_root_dir + var_name + "/" + data_name + '_' + pat + '_scapny.h5ad'):
        adata = sc.read_h5ad(data_root_dir + var_name + "/" + data_name + '_' + pat + '_scapny.h5ad')
    else:
        return "apple"

    ContrastiveVI.setup_anndata(adata, layer="counts")

    if perturb_only == "yes":
        target_adata = adata[adata.obs["expand"] == "yes"].copy()
        background_adata = adata[adata.obs["expand"] == "no"].copy()

        model = ContrastiveVI(
            adata,
            n_salient_latent=10,
            n_background_latent=10,
            use_observed_lib_size=False
        )
    
        background_indices = np.where(adata.obs["expand"] == "no")[0]
        target_indices = np.where(adata.obs["expand"] == "yes")[0]
    else:
        target_adata = adata[adata.obs["condition"] != "unperturb"].copy()
        background_adata = adata[adata.obs["condition"] == "unperturb"].copy()

        model = ContrastiveVI(
            adata,
            n_salient_latent=10,
            n_background_latent=10,
            use_observed_lib_size=False
        )

        background_indices = np.where(adata.obs["condition"] == "unperturb")[0]
        target_indices = np.where(adata.obs["condition"] != "unperturb")[0]


    model.train(
        check_val_every_n_epoch=1,
        train_size=0.8,
        background_indices=background_indices,
        target_indices=target_indices,
        use_gpu=True,
        early_stopping=True,
        max_epochs=500,
    )

    target_adata.obsm['salient_rep'] = model.get_latent_representation(target_adata, representation_kind='salient')
    adata.obsm['background_rep'] = model.get_latent_representation(adata, representation_kind='background')
    adata.obsm['salient_rep'] = model.get_latent_representation(adata, representation_kind='salient')
    
    
    adata.write(data_name + "_" + pat + "_contrast.h5ad")
    #target_adata = adata[adata.obs["Tissue"] != "Normal"].copy()
    #target_adata.obsm['salient_rep'] = model.get_latent_representation(target_adata, representation_kind='salient')
    
    
    k = 5
    if adata.obs["expand"].nunique() == 1:
        return "apple"

    if perturb_only == "yes":
        latent_basal_entropy  = entropy_process(adata, "background_rep", 'expand')
        latent_after_entropy = entropy_process(adata, "salient_rep", 'expand')
        expand_entropy = entropy_process(adata, "salient_rep","expand")

        statistic, p_value = wilcoxon(latent_basal_entropy, latent_after_entropy,alternative='two-sided')

        if expand_only == "yes":
            target_adata1 = target_adata[target_adata.obs["expand"] == "yes"]
        else:
            target_adata1 = target_adata.copy()


    else:
        latent_basal_entropy  = entropy_process(adata, "background_rep", 'condition')
        latent_after_entropy = entropy_process(adata, "salient_rep", 'condition')
        expand_entropy = entropy_process(target_adata, "salient_rep","expand")

        statistic, p_value = wilcoxon(latent_basal_entropy, latent_after_entropy,alternative='two-sided')

        if expand_only == "yes":
            target_adata1 = target_adata[(target_adata.obs["expand"].values == "yes") & (target_adata.obs["condition"].values == "perturb")]
        else:
            target_adata1 = target_adata[target_adata.obs["condition"] == "perturb"]

    if target_adata1.n_obs < 50:
        corr_res = ["NA", "NA"]
    else:
        corr_res = clonality_correlation(target_adata1, k)

    
    results = [latent_basal_entropy.mean(), latent_after_entropy.mean(), p_value, expand_entropy.mean(), corr_res[0], corr_res[1]]
    return(results)


# In[6]:

pat_list = list(line.strip() for line in open(sample_list_file))

if not os.path.exists(save_root_dir + var_name):
    os.makedirs(save_root_dir + var_name)
os.chdir(save_root_dir + var_name)
total_result = pd.DataFrame()

for pat in pat_list:
    print(pat)
    pat_res = contrastive_process(pat)
    if pat_res == "apple":
        continue

    df = pd.DataFrame(pat_res)
    df.index = ["latent_basal_entropy","latent_after_entropy","p_value","expand_entropy","expand_corr", "corr_pval"]
    total_result[pat] = pd.Series(pat_res)
total_result.index = ["latent_basal_entropy","latent_after_entropy","p_value","expand_entropy","expand_corr", "corr_pval"]
total_result.to_csv(data_name + "_contrastive_result", sep = "\t")


# In[ ]:




