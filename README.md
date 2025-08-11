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
res_plot(adata, "cpa_perturb_only", latent = None)<br/><br/><br/>
