# BRiCE (Benchmarking RespondIng cells by Clonal Expansion)
Scripts and tutorial for BRiCE<br/>

Relevant paper: Benchmarking of predicting the responding cells upon perturbation affecting both gene expression and cellular abundance in scRNA-sequencing <br/>

Citation: not yet<br/>

##
Folder description<br/>
BRiCE: contains functions for BRiCE<br/>
example_data: example_data<br/>
manuscrip_methods: scripts used for the manuscript<br/>

##
Tutorial<br/><br/>

Preprocess<br/>
R
#library(Seurat)
#library(SingleCellExperiment)
#library(zellkonverter)

source("preprocess.R")<br/>


#meta information<br/>
#orig.ident<br/>
#clonotype<br/>
#individual<br/>
#condition: unperturb, perturb<br/>

clonal_seurat_process(user_seurat_object, "user_file_name")<br/>
-This will generate "h5ad" file <br/>
##
BRiCE<br/>
python: 3.9<br/>

import sys, os<br/>
import scanpy as sc<br/>
import pandas as pd<br/>
import numpy as np<br/>
import scipy.sparse as sp<br/>
from scipy.stats import spearmanr<br/>
import argparse<br/>
from math import log<br/>
import random<br/>
import bisect<br/>
from scipy.stats import wilcoxon<br/>
import seaborn as sns<br/>
import matplotlib.pyplot as plt<br/>
import glob<br/>

sys.path.append("/path/to/performance.py")<br/>
sys.path.append("/path/to/visualization.py")<br/>


adata = sc.read_h5ad("bcr_vaccine_hvg_only.h5ad")<br/>
obj1 = adata[(adata.obs["condition"] == "perturb")]<br/>
result = entropy_process(obj1, "expand", latent = None)<br/>
result.mean()<br/>
#0.05934534403705289<br/>

obj = adata[(adata.obs["condition"] == "perturb") & (adata.obs["expand"] == "yes")]<br/>
result = clonality_correlation(obj, k=5, latent = None)<br/>
result[0]<br/>
#0.06120139489080349<br/>

rc_umap(adata, "test", "X_umap", need_neighbor = False)<br/>
res_plot(adata, "test", latent = None)<br/>

<img src="/example_data/test_RC_umap.png" width="600"/>
<img src="/example_data/test_Res_plot.png" width="600"/>

##
adata = sc.read_h5ad("bcr_vaccine_vaccin_cpa_after.h5ad")<br/>
obj1 = adata[(adata.obs["condition"] == "perturb")]<br/>
result = entropy_process(obj1, "expand", latent = None)<br/>
result.mean()<br/>
#0.0<br/>

obj = adata[(adata.obs["condition"] == "perturb") & (adata.obs["expand"] == "yes")]<br/>
result = clonality_correlation(obj, k=5, latent = None)<br/>
result[0]<br/>
#0.018688870657997153<br/>

rc_umap(adata, "cpa_perturb_only", "X_umap", need_neighbor = True)<br/>
res_plot(adata, "cpa_perturb_only", latent = None)<br/>

<img src="/example_data/cpa_perturb_only_RC_umap.png" width="600"/>
<img src="/example_data/cpa_perturb_only_Res_plot.png" width="600"/>
