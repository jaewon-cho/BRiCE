#!/usr/bin/env python
# coding: utf-8

# In[1]:
#usage: scgen_run.py sample_list_file data_name data_root_dir var_name save_root_dir
#sample_file_list:
#data_name: gse236581
# data_root_dir include last  /
# /home/dreadcupper/jaewon/perturb_clonal_method/anaylsis/gse236581/common_python/
#var_name: var2000
# save_root_dir include last  /
# /home/dreadcupper/jaewon/perturb_clonal_method/anaylsis/gse236581/scgen/


import sys, os
import logging
import scanpy as sc
import scgen

import pandas as pd

import argparse
from math import log
import random
import bisect
from scipy.stats import wilcoxon
#from scvi.data import setup_anndata

sample_list_file = sys.argv[1]
data_name = sys.argv[2]
data_root_dir = sys.argv[3]
var_name = sys.argv[4]
save_root_dir = sys.argv[5]

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


# In[4]:


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


# In[5]:


def scgen_process(pat):
    if os.path.isfile(data_root_dir + var_name + "/" + data_name + '_' + pat + '_scapny.h5ad'):
        adata = sc.read_h5ad(data_root_dir + var_name + "/" + data_name + '_' + pat + '_scapny.h5ad')
    else:
        return "apple"

    scgen.SCGEN.setup_anndata(adata, batch_key="condition")
    model = scgen.SCGEN(adata)
    model.train(
        max_epochs=100,
        batch_size=32,
        early_stopping=True,
        early_stopping_patience=25
    )
    latent_X = model.get_latent_representation()
    latent_adata = sc.AnnData(X=latent_X, obs=adata.obs.copy())
    latent_adata.write(data_name + "_" + pat + "_scgen.h5ad")
    
    latent_after_entropy = entropy_process(latent_adata,'condition')
    
    new_obj1 = latent_adata[latent_adata.obs["condition"] == "perturb"]
    expand_entropy = entropy_process(new_obj1, "expand")
    
    #statistic, p_value = wilcoxon(latent_basal_entropy, latent_after_entropy,alternative='two-sided')
    sc.pp.neighbors(new_obj1, n_neighbors=30, knn = True)
    random_result = []
    for i in range(1000):
        #print(i)
        random_exp_entropy = entropy_rand_process(new_obj1, "expand")
        tmp_val = random_exp_entropy.mean()
        random_result.append(tmp_val)
    
    random_result=sorted(random_result)
    random_mean = sum(random_result)/len(random_result)
    expand_rank = bisect.bisect_left(random_result, expand_entropy.mean()) + 1
    #empirical p value less
    
    
    results = [latent_after_entropy.mean(), expand_entropy.mean(), random_mean, expand_rank/1000]
    return(results)


# In[9]:

pat_list = list(line.strip() for line in open(sample_list_file))

if not os.path.exists(save_root_dir + var_name):
    os.makedirs(save_root_dir + var_name)
os.chdir(save_root_dir + var_name)
total_result = pd.DataFrame()

for pat in pat_list:
    print(pat)
    pat_res = scgen_process(pat)
    if pat_res == "apple":
        continue
    #sinkhorn_knopp does not fit
    
    df = pd.DataFrame(pat_res)
    df.index = ["latent_after_entropy","expand_entropy","random_entropy","exp_pval"]
    total_result[pat] = pd.Series(pat_res)
total_result.index = ["latent_after_entropy","expand_entropy","random_entropy","exp_pval"]
total_result.to_csv(data_name + "_scgen_result", sep = "\t")


# In[ ]:




