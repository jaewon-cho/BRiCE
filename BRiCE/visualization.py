import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

#adata: object
#name: output file name
#embed: umap embedding key: scanpy: "X_umap"
#need_neighbor: if sc.pp.neighbor needed: True (default: False)

def rc_umap(adata, name, embed, need_neighbor = False):
 #adata = sc.read_h5ad("/home/dreadcupper/jaewon/perturb_clonal_method/anaylsis/10x_bcell_vaccin/common_python/var2000/bcr_vaccine_vaccin_scapny.h5ad")
 if need_neighbor:
      sc.pp.neighbors(adata, knn = True)
      sc.tl.umap(adata)
      embed = "X_umap"
 ref_tumor = adata[adata.obs['condition'] == 'perturb'].copy()
 ref_tumor.obs['expand'] = pd.Categorical(ref_tumor.obs['expand'], categories=["no", "yes"], ordered=True)
 ref_tumor = ref_tumor[ref_tumor.obs.sort_values("expand").index]
# Extract UMAP embedding (assuming it’s already computed)
 umap = ref_tumor.obsm[embed]  # or 'embed' if that’s your key
# Create plot DataFrame
 plot_df = pd.DataFrame({
    "UMAP1": umap[:, 0],
    "UMAP2": umap[:, 1],
    "Expand": ref_tumor.obs['expand'].values
 })
# Set factor order and colors
 plot_df['Expand'] = pd.Categorical(plot_df['Expand'], categories=["yes", "no"], ordered=True)
 colors = {"yes": "#F8B62D", "no": "#595757"}
# Plot
 plt.figure(figsize=(6, 6), dpi=450)
 sns.scatterplot(
    data=plot_df,
    x="UMAP1", y="UMAP2",
    hue="Expand",
    palette=colors,
    s=30, edgecolor=None
 )
 plt.xlabel("UMAP1")
 plt.ylabel("UMAP2")
 sns.despine()
 plt.legend(title="Expand", loc="upper right")
 plt.tight_layout()
# Save to file
 plt.savefig(f"{name}_RC_umap.png")
 plt.close()

#adata: object
#name: output file name
#latent: use_rep for sc.pp.neighbors
#k = number of neighbors for kNN

def res_plot(adata, name, latent = None, k =5):
 #adata = sc.read_h5ad("/home/dreadcupper/jaewon/perturb_clonal_method/anaylsis/gse236581/common_python/var2000/gse236581_P20_scapny.h5ad")
 ref_tumor = adata[(adata.obs['condition'] == 'perturb') & (adata.obs['expand'] == 'yes')].copy()
 ref_tumor.obs['expand'] = pd.Categorical(ref_tumor.obs['expand'], categories=["no", "yes"], ordered=True)
 ref_tumor = ref_tumor[ref_tumor.obs.sort_values("expand").index]
 sc.pp.neighbors(ref_tumor, use_rep = latent, n_neighbors=k+1, knn = True)
 knn_binary = (ref_tumor.obsp["distances"] > 0).astype(int)
 clonality = ref_tumor.obs['clonality'].values
 neighbor_sums = knn_binary@clonality
 neighbor_counts = np.array(knn_binary.sum(1)).flatten()
 mean_clonality = neighbor_sums / neighbor_counts
# Plot
 x = clonality
 y = mean_clonality
 x = np.log2(clonality)
 y = np.log2(mean_clonality)
 plt.figure(figsize=(6, 6), dpi=450)
 plt.scatter(x, y, s=15)
 plt.xlabel("Log2(Clonality)")
 plt.ylabel("Log2(Mean Neighbor Clonality)")
 plt.tight_layout()
# Save to file
 plt.savefig(f"{name}_Res_plot.png")
 plt.close()


