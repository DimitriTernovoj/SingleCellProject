import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb

#read in anndata object & annotations
adata = sc.read_h5ad(snakemake.input[0])
annotation_list = np.load(snakemake.input[1])
annotation_list = annotation_list.tolist()

if snakemake.params.cluster_method == "specter":
  current = "specter"
else:
  current = "louvain_r"+ snakemake.params.clustering_resolution


#exchange "unannotateds" with names specified by the user (if selected) + visualization - names must be as long as there are unannotateds in data
names = {snakemake.params.names_for_unannotated}.pop().split(",")
if len(names[0]) != 0:
  count = 0
  for i in range(len(annotation_list)):
    if annotation_list[i][0:11] == "unannotated":
       annotation_list[i] = names[count]
       count += 1

  adata.rename_categories(current, annotation_list)
  sc.pl.umap(adata, color=current, size=15, legend_loc='on data', save="_named_annotated_clustering")


#prepare subclustering defined by user
further_subclusterings = {snakemake.params.further_subclusterings}.pop().split(",")
temp = []

for i in range(len(annotation_list)):
  for j in range(len(further_subclusterings)):
    if annotation_list[i] == further_subclusterings[j]:
      temp.append(True)
      break
  if len(temp)!= i+1:
    temp.append(False)

#subclustering
for i in range(len(annotation_list)):
  if len(annotation_list[i].split("_")) != 1 or temp[i]:
    sc.tl.louvain(adata, restrict_to=(current, [annotation_list[i]]), resolution=float({snakemake.params.subclustering_resolution}.pop()), key_added=current+str(i))
    current=current+str(i)

if current+"_colors" in adata.uns:
  del adata.uns[current+"_colors"]

#visualization
sc.pl.umap(adata, color=current, palette=sc.pl.palettes.godsnot_64, save="_final_clustering")
sc.pl.umap(adata, color='region', palette=sc.pl.palettes.godsnot_64, save="_regions_on_finalclustering")


#what clusters to include for further analysis
clusters = snakemake.params.clusters_to_include.split(",")
clusters_to_include = []

if len(clusters) != 0:
  for i in adata.obs[current].cat.categories:
    for j in clusters:
      if i.startswith(j):
        clusters_to_include.append(i)
else:
  clusters_to_include =  adata.obs[current].cat.categories


#prepare adata_ent object for slingshot (trajectory inference)
adata_ent = adata[np.isin(adata.obs[current], clusters_to_include),:].copy()

#subset to highly variable genes
sc.pp.highly_variable_genes(adata_ent, flavor='seurat', n_top_genes= int(snakemake.params.amount_of_hvgs), subset=False)

#recalculating PCA for subset
sc.pp.pca(adata_ent, svd_solver='arpack', use_highly_variable=True, n_comps=7)
sc.pl.pca(adata_ent, save="_subset_for_ti")
sc.pl.pca_variance_ratio(adata_ent, save="_for_ti")

adata_ent.write("file_dir/adata_ent.h5ad")

#plot diffusion pseudotime on adata_ent
sc.pp.neighbors(adata_ent)
sc.tl.diffmap(adata_ent)

sc.pl.diffmap(adata_ent, components='1,2', color=current, save="_diffusion_pseudotime_1")
sc.pl.diffmap(adata_ent, components='1,3', color=current, save="_diffusion_pseudotime_2")


#prepare adata_ent_nbc object for slingshot (trajectory inference) with non-batch-corrected data
#subsetting data set - non-batch corrected
cell_mask = np.isin(adata.obs[current], clusters_to_include)
adata_ent_nbc = sc.AnnData(X = adata.raw.X[cell_mask,:])
adata_ent_nbc.obs = adata.obs[cell_mask]
adata_ent_nbc.var = adata.var.copy()

#subset to highly variable genes
sc.pp.highly_variable_genes(adata_ent_nbc, flavor='cell_ranger', n_top_genes=int(snakemake.params.amount_of_hvgs), subset=True)

#recalculating PCA for subset
sc.pp.pca(adata_ent_nbc, svd_solver='arpack')
sc.pl.pca(adata_ent_nbc, save="_subset_for_ti_nbc")
sc.pl.pca_variance_ratio(adata_ent_nbc, save="_for_ti_nbc")

adata_ent_nbc.obsm['X_pca'] = adata_ent_nbc.obsm['X_pca'][:,0:7] #hier die Zahl variabel machen!
adata_ent_nbc.X = sp.sparse.csr_matrix.todense(adata_ent_nbc.X)

#export as separate files
data_mat3 = adata_ent_nbc.X.T
np.savetxt("file_dir/data_mat3.csv", data_mat3, delimiter=",")
adata_ent_nbc.obs.to_csv("file_dir/obs_nbc.csv")
adata_ent_nbc.var.to_csv("file_dir/var_nbc.csv")
np.savetxt("file_dir/pca_nbc.csv", adata_ent_nbc.obsm["X_pca"], delimiter=",")


#preperation for monocle2
data_mat_mon = adata.layers['counts'].T
var_mon=adata.var.copy()
obs_mon=adata.obs.copy()

np.savetxt("file_dir/data_mat_mon.csv", data_mat_mon, delimiter=",")
obs_mon.to_csv("file_dir/obs_mon.csv")
var_mon.to_csv("file_dir/var_mon.csv")


#prepare for differential testing
adata_dt = adata.copy()
adata_dt.X = adata.raw.X
adata_dt.obs['n_genes'] = (adata_dt.X > 0).sum(1)
adata_dt.X = sp.sparse.csr_matrix(adata_dt.X)

#choose correct file based on whether sphetcher was used or not
if snakemake.params.normed_type == "sphetcher":
  name = "file_dir/normed_data_sphetcher.h5ad"
else:
  name= "file_dir/normed_data_normal.h5ad"

normed_adata = sc.read_h5ad(name)
normed_adata.X = sp.sparse.csr_matrix(normed_adata.X)
adata_dt.raw = normed_adata

adata_dt.write("file_dir/adata_dt.h5ad")


#safe final version of anndata object
adata.write("file_dir/post_clustering_final.h5ad")
