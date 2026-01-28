#!/usr/bin/env python
# coding: utf-8

# In[13]:
#usage: cpa_run.py sample_list_file data_name data_root_dir var_name save_root_dir perturb_only expand_only
#sample_file_list:
#data_name: gse236581
# data_root_dir include last  /
# /home/dreadcupper/jaewon/perturb_clonal_method/anaylsis/gse236581/common_python/
#var_name: var2000
# save_root_dir include last  /
# /home/dreadcupper/jaewon/perturb_clonal_method/anaylsis/gse236581/cpa/
#perturb_only: yes or no
#yes: perturb only
#expand_only: yes or no
#yes: expand only for measuring clonality correlation



#Predicting perturbation responses for unseen cell-types (context transfer)
import sys, os
import cpa
import scanpy as sc
import pandas as pd
import numpy as np

import argparse
from math import log
import random
import bisect
from scipy.stats import wilcoxon
from scipy.stats import spearmanr

# In[2]:
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


# In[3]:


def entropy_process(latent_obj, category):
    #category: expand or Tissue
    dataset=latent_obj
    sc.pp.neighbors(dataset, n_neighbors=31, knn = True)
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
    sc.pp.neighbors(adata1, use_rep = latent, n_neighbors=k+1, knn = True)
    knn_binary = (adata1.obsp["distances"] > 0).astype(int)
    clonality = adata1.obs['clonality'].values
    neighbor_sums = knn_binary@clonality
    neighbor_counts = np.array(knn_binary.sum(1)).flatten()
    mean_clonality = neighbor_sums / neighbor_counts

    corr, p_val = spearmanr(clonality, mean_clonality)
    return([corr, p_val])


# In[12]:


def CPA_process(pat):
    if os.path.isfile(data_root_dir + var_name + "/" + data_name + '_' + pat + '_scapny.h5ad'):
        adata = sc.read_h5ad(data_root_dir + var_name + "/" + data_name + '_' + pat + '_scapny.h5ad')
    else:
        return "apple"

    if perturb_only == "yes":
        cpa.CPA.setup_anndata(adata,
                          perturbation_key='expand',
                          control_group='no',
                          is_count_data=False,
                          max_comb_len=1,
                         )
    else:
        cpa.CPA.setup_anndata(adata,
                          perturbation_key='condition',
                          control_group='unperturb',
                          is_count_data=False,
                          max_comb_len=1,
                         )

    
    model = cpa.CPA(adata=adata,
                **model_params,
               )
    
    model.train(max_epochs=2000,
            use_gpu=True,
            batch_size=512,
            plan_kwargs=trainer_params,
            early_stopping_patience=5,
            check_val_every_n_epoch=5,
            save_path='cpa_save',
           )
    
    latent_outputs = model.get_latent_representation(adata, batch_size=2048)
    basal_obj = latent_outputs['latent_basal']
    after_obj = latent_outputs['latent_after']
    basal_obj.write(data_name + "_" + pat + "_cpa_basal.h5ad")
    after_obj.write(data_name + "_" + pat + "_cpa_after.h5ad")

    k = 5
    
    if perturb_only == "yes":
        latent_basal_entropy = entropy_process(latent_outputs['latent_basal'], "expand")
        latent_after_entropy = entropy_process(latent_outputs['latent_after'], "expand")
        new_obj = latent_outputs['latent_after']
        new_obj1 = new_obj.copy()
        expand_entropy = entropy_process(new_obj1, "expand")

        statistic, p_value = wilcoxon(latent_basal_entropy, latent_after_entropy,alternative='two-sided')
        if expand_only == "yes":
            new_obj11 = new_obj1[new_obj1.obs["expand"] == "yes"]
        else:
            new_obj11 = new_obj1.copy()


    else:
        latent_basal_entropy = entropy_process(latent_outputs['latent_basal'], "condition")
        latent_after_entropy = entropy_process(latent_outputs['latent_after'], "condition")
        new_obj = latent_outputs['latent_after']
        new_obj1 = new_obj[new_obj.obs["condition"] == "perturb"]
        expand_entropy = entropy_process(new_obj1, "expand")

        statistic, p_value = wilcoxon(latent_basal_entropy, latent_after_entropy,alternative='two-sided')
        if expand_only == "yes":
            new_obj11 = new_obj1[(new_obj1.obs["expand"].values == "yes") & (new_obj1.obs["condition"].values == "perturb")]
        else:
            new_obj11 = new_obj1[new_obj1.obs["condition"] == "perturb"]

    if new_obj11.n_obs < 50:
        corr_res = ["NA", "NA"]
    else:
        corr_res = clonality_correlation(new_obj11, k)
    
    
    results = [latent_basal_entropy.mean(), latent_after_entropy.mean(), p_value, expand_entropy.mean(), corr_res[0], corr_res[1]]
    return(results)


# In[6]:


model_params = {
    "n_latent": 64,
    "recon_loss": "nb",
    "doser_type": "linear",
    "n_hidden_encoder": 128,
    "n_layers_encoder": 2,
    "n_hidden_decoder": 512,
    "n_layers_decoder": 2,
    "use_batch_norm_encoder": True,
    "use_layer_norm_encoder": False,
    "use_batch_norm_decoder": False,
    "use_layer_norm_decoder": True,
    "dropout_rate_encoder": 0.0,
    "dropout_rate_decoder": 0.1,
    "variational": False,
    "seed": 6977,
}

trainer_params = {
    "n_epochs_kl_warmup": None,
    "n_epochs_pretrain_ae": 30,
    "n_epochs_adv_warmup": 50,
    "n_epochs_mixup_warmup": 0,
    "mixup_alpha": 0.0,
    "adv_steps": None,
    "n_hidden_adv": 64,
    "n_layers_adv": 3,
    "use_batch_norm_adv": True,
    "use_layer_norm_adv": False,
    "dropout_rate_adv": 0.3,
    "reg_adv": 20.0,
    "pen_adv": 5.0,
    "lr": 0.0003,
    "wd": 4e-07,
    "adv_lr": 0.0003,
    "adv_wd": 4e-07,
    "adv_loss": "cce",
    "doser_lr": 0.0003,
    "doser_wd": 4e-07,
    "do_clip_grad": True,
    "gradient_clip_value": 1.0,
    "step_size_lr": 10,
}


# In[8]:

pat_list = list(line.strip() for line in open(sample_list_file))

if not os.path.exists(save_root_dir + var_name):
    os.makedirs(save_root_dir + var_name)
os.chdir(save_root_dir + var_name)
total_result = pd.DataFrame()

for pat in pat_list:
    print(pat)
    pat_res = CPA_process(pat)
    if pat_res == "apple":
        continue

    df = pd.DataFrame(pat_res)
    df.index = ["latent_basal_entropy","latent_after_entropy","p_value","expand_entropy","expand_corr", "corr_pval"]
    total_result[pat] = pd.Series(pat_res)
total_result.index = ["latent_basal_entropy","latent_after_entropy","p_value","expand_entropy","expand_corr", "corr_pval"]
total_result.to_csv(data_name + "_cpa_result", sep = "\t")


# In[ ]:




