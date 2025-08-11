#!/usr/bin/env python
# coding: utf-8

# In[2]:

#usage: ./cinema_ot_run.py sample_list_file data_name data_root_dir var_name save_root_dir

#sample_file_list: 
#data_name: gse236581
# data_root_dir include last  /
# /home/dreadcupper/jaewon/perturb_clonal_method/anaylsis/gse236581/common_python/
#var_name: var2000
# save_root_dir include last  /
# /home/dreadcupper/jaewon/perturb_clonal_method/anaylsis/gse236581/cinema_ot/

import numpy as np
import scanpy as sc
import cinemaot as co
import matplotlib.colors as colors
import matplotlib.pyplot as plt

import sys,os
sys.path.append("/home/dreadcupper/.conda/envs/cinema_ot/lib/python3.9/site-packages/cinemaot")
import utils
import pandas as pd

import argparse
from math import log
import random
import bisect
from scipy.stats import wilcoxon

sample_list_file = sys.argv[1]
data_name = sys.argv[2]
data_root_dir = sys.argv[3]
var_name = sys.argv[4]
save_root_dir = sys.argv[5]


# In[3]:


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


# In[4]:


def entropy_process(latent_obj, latent_category, category):
    #category: expand or Tissue
    dataset=latent_obj
    sc.pp.neighbors(dataset, use_rep= latent_category, n_neighbors=30, knn = True)
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


# In[5]:


def entropy_rand_process(dataset, category):
    #dataset=latent_obj
    #
    #n_neighbor: batchbench
    
    knn_graph = dataset.obsp['connectivities']
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


# In[6]:


def cinema_process(pat):
    if os.path.isfile(data_root_dir + var_name + "/" + data_name + '_' + pat + '_scapny.h5ad'):
        adata = sc.read_h5ad(data_root_dir + var_name + "/" + data_name + '_' + pat + '_scapny.h5ad')
    else:
        return "apple"
    
    cf, ot, de = co.cinemaot.cinemaot_unweighted(adata,obs_label='condition', ref_label='perturb', expr_label='unperturb',mode='parametric',thres=0.5,smoothness=1e-5,eps=1e-3)
    if True in np.isnan(ot):
        return "apple"
    
    adata.obsm['cf'] = cf.copy()
    adata.obsm['cf'][adata.obs['condition']=='perturb',:] = np.matmul(ot/np.sum(ot,axis=1)[:,None],cf[adata.obs['condition']=='unperturb',:])
    sc.pp.pca(de)
    
    adata.write(data_name + "_" + pat + "_cinema_cf.h5ad")
    de.write(data_name + "_" + pat + "_cinema_de.h5ad")
    np.save(data_name + "_" + pat + "_cinema_ot.h5ad", ot)
    #target_adata = adata[adata.obs["Tissue"] != "Normal"].copy()
    #target_adata.obsm['salient_rep'] = model.get_latent_representation(target_adata, representation_kind='salient')
    
    latent_basal_entropy  = entropy_process(adata, "cf", 'condition')
    #latent_after_entropy = entropy_process(adata, "salient_rep", 'Tissue')
    expand_entropy = entropy_process(de, "X_embedding","expand")
    
    #statistic, p_value = wilcoxon(latent_basal_entropy, latent_after_entropy,alternative='two-sided')
    sc.pp.neighbors(de, use_rep= "X_embedding", n_neighbors=30, knn = True)
    random_result = []
    for i in range(1000):
        #print(i)
        random_exp_entropy = entropy_rand_process(de, "expand")
        tmp_val = random_exp_entropy.mean()
        random_result.append(tmp_val)
    
    random_result=sorted(random_result)
    random_mean = sum(random_result)/len(random_result)
    expand_rank = bisect.bisect_left(random_result, expand_entropy.mean()) + 1
    #empirical p value less
    
    
    results = [latent_basal_entropy.mean(), expand_entropy.mean(), random_mean, expand_rank/1000]
    return(results)


# In[ ]:


pat_list = list(line.strip() for line in open(sample_list_file))
total_result = pd.DataFrame()

if not os.path.exists(save_root_dir + var_name):
    os.makedirs(save_root_dir + var_name)
os.chdir(save_root_dir + var_name)

for pat in pat_list:
    print(pat)
    pat_res = cinema_process(pat)
    if pat_res == "apple":
        continue
    #sinkhorn_knopp does not fit
    
    df = pd.DataFrame(pat_res)
    df.index = ["latent_basal_entropy","expand_entropy","random_entropy","exp_pval"]
    total_result[pat] = pd.Series(pat_res)
total_result.index = ["latent_basal_entropy","expand_entropy","random_entropy","exp_pval"]
total_result.to_csv(data_name + "_cinema_result", sep = "\t")


# cf: original confounder embedding

# de: effective confounder embedding

# In[ ]:





# In[ ]:





# In[ ]:




