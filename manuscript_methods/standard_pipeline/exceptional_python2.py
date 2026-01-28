#!/usr/bin/env python
# coding: utf-8
#usage: ./exceptional_python.py data_name data_name_extra data_root_dir var_name 
#data_name: gse236581
#data_name_extra = tcell_seurat
#data_dir:/home/dreadcupper/jaewon/perturb_clonal_method/anaylsis/gse236581/
#var_name: perturb_deg, expand_deg, corr 


import sys, os
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

data_name = sys.argv[1]
data_name_extra = sys.argv[2]
data_dir = sys.argv[3]
var_name = sys.argv[4]


adata = sc.read_h5ad(data_dir + data_name + "_" + data_name_extra +  '.h5ad')


if not os.path.exists(data_dir + "common_python/" + var_name):
    os.makedirs(data_dir + "common_python/" + var_name)
os.chdir(data_dir + "common_python/" + var_name)

pat_list = set(adata.obs["individual"])


for pat in pat_list:
    print(pat)
    adata1 = adata[adata.obs["individual"] == pat]
    adata1.layers["counts"] = adata1.X.copy()
    sc.pp.normalize_total(adata1)
    sc.pp.log1p(adata1)
    sc.pp.highly_variable_genes(adata1, n_top_genes=500)
    if var_name == "perturb_deg":
        sc.tl.rank_genes_groups(adata1, groupby='condition', method='wilcoxon')
        deg_df = sc.get.rank_genes_groups_df(adata1, group='perturb')
        sig_genes =  deg_df.loc[deg_df['pvals_adj'] < 0.05, 'names'].tolist()
        new_var_gn = set(sig_genes)
    elif var_name == "expand_deg":
        adata2 = adata1[adata1.obs['condition'] == 'perturb'].copy()
        sc.tl.rank_genes_groups(adata2, groupby='expand', method='wilcoxon')
        deg_df = sc.get.rank_genes_groups_df(adata2, group='yes')
        valid_pvals = deg_df['pvals']
        _, qvals, _, _ = multipletests(valid_pvals, method='fdr_bh')
        deg_df['qval'] = qvals
        sig_genes =  deg_df.loc[deg_df['qval'] < 0.1, 'names'].tolist()
        new_var_gn = set(sig_genes)
    elif var_name == "corr":
        adata2 = adata1[adata1.obs['condition'] == 'perturb'].copy()
        X = adata2.X
        if not isinstance(X, np.ndarray):
            X = X.toarray()

        clonality = adata2.obs['clonality']
        gene_names = adata2.var_names

        # Create a DataFrame with gene expression and clonality
        expr_df = pd.DataFrame(X, columns=gene_names)
        expr_df['clonality'] = clonality.values

    # Group by clonality and average gene expression
        grouped = expr_df.groupby('clonality').mean()

# Remove 'clonality' if still in columns after groupby().mean()
        if 'clonality' in grouped.columns:
            grouped = grouped.drop(columns='clonality')

# Compute Spearman correlation
        clonality_values = grouped.index.values
        correlations = {}
        pvals = {}

        for gene in grouped.columns:
            rho, pval = spearmanr(clonality_values, grouped[gene].values)
            correlations[gene] = rho
            pvals[gene] = pval

    # Compile results
        results_df = pd.DataFrame({
            'spearman_rho': correlations,
            'pval': pvals
        })
        
        results_df['qval'] = np.nan
        valid_mask = results_df['pval'].notna()
        valid_pvals = results_df.loc[valid_mask, 'pval']

        _, qvals, _, _ = multipletests(valid_pvals, method='fdr_bh')
        results_df.loc[valid_mask, 'qval'] = qvals
        sig_db = results_df[results_df['qval'] < 0.1].sort_values('pval')
        new_var_gn = sig_db.index.tolist()

    if len(new_var_gn) <= 30:
        continue

    adata1.var['highly_variable'] = adata1.var_names.isin(new_var_gn)
    sc.tl.pca(adata1)
    sc.pp.neighbors(adata1, n_neighbors=10, n_pcs=30)
    sc.tl.umap(adata1)
    adata1.write(data_name + "_" + pat + "_scapny.h5ad")


