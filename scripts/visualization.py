import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb

#read in anndata object
adata = sc.read_h5ad(snakemake.input[0])

#calculate pca + distances between the cells
sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adata)

#calculate visualizations
sc.tl.tsne(adata, n_jobs=12)
sc.tl.umap(adata)
sc.tl.diffmap(adata)
sc.tl.draw_graph(adata)

#plot the visualizations
sc.pl.pca_scatter(adata, color='n_counts', save="_datavis_" + snakemake.wildcards.downsampling)
sc.pl.tsne(adata, color='n_counts', save="_datavis_"+ snakemake.wildcards.downsampling)
sc.pl.umap(adata, color='n_counts', save="_datavis_"+ snakemake.wildcards.downsampling)
sc.pl.diffmap(adata, color='n_counts', save="_datavis_" + snakemake.wildcards.downsampling)
sc.pl.draw_graph(adata, color='n_counts', save="_datavis_"+ snakemake.wildcards.downsampling)

#visualizing cell cycle effects
if snakemake.params.ref_genes != "":
  #read in cell cycle genes
  cc_genes = pd.read_table({snakemake.params.ref_genes}.pop())

  s_genes_mm = cc_genes['s'].dropna()
  g2m_genes_mm = cc_genes['g2m'].dropna()

  s_genes_mm_ens = adata.var_names[np.in1d(adata.var_names, s_genes_mm)]
  g2m_genes_mm_ens = adata.var_names[np.in1d(adata.var_names, g2m_genes_mm)]

  #calculate cell cycle scoring
  sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes_mm_ens, g2m_genes=g2m_genes_mm_ens)

  #plot the visualizations
  sc.pl.umap(adata, color=['S_score', 'G2M_score'], use_raw=False, save="_cellcycle_1_" + snakemake.wildcards.downsampling)
  sc.pl.umap(adata, color='phase', use_raw=False, save="_cellcycle_2_" + snakemake.wildcards.downsampling)

  #possible correcting of cell cycle effects (not tested yet)
  #sc.pp.regress_out(adata,["S_score","G2M_score"])
  #sc.pp.scale(adata)

name = "file_dir/post_vis_data_" + snakemake.wildcards.downsampling + ".h5ad"
adata.write(name)

#saving files for clustering with specter
#np.savetxt("file_dir/matlab_test_X.csv", adata.X, delimiter=",")
if snakemake.params.clustering == "specter":
  np.savetxt("file_dir/specter_pca_data_" + snakemake.wildcards.downsampling + ".csv", adata.obsm["X_pca"], delimiter=",")

