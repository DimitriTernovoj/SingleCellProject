import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb

#read in specter results and the anndata object
specter_clustering = pd.read_csv("file_dir/specter_clustering.csv", sep=",", header=None)
specter_clustering = specter_clustering.to_numpy()
adata = sc.read_h5ad(snakemake.input[0])

#perform clustering - using highly variable genes
sc.tl.louvain(adata, key_added='louvain_r1')

#add specter-clustering to anndata.obs
adata.obs["specter"] = specter_clustering
adata.obs["specter"] = adata.obs["specter"].astype("category")

#visualize the clustering and how this is reflected by different technical covariates
sc.pl.umap(adata, color=['louvain_r1', 'specter'], palette=sc.pl.palettes.vega_20, save="_clustering_result")
sc.pl.umap(adata, color=['region', 'n_counts'], save="_region_ncounts_oncluster")
sc.pl.umap(adata, color=['log_counts', 'mt_frac'], save="_logcounts_mtfrac_oncluster")

#marker genes & cluster annotation
#calculate marker genes
sc.tl.rank_genes_groups(adata, groupby='specter', key_added='rank_genes_specter')

#Plot marker genes
sc.pl.rank_genes_groups(adata, key='rank_genes_specter', fontsize=12, save="_markergenes")

#read in specified marker genes
temp = pd.read_table({snakemake.params.celltype_annotation}.pop())
marker_genes = {}
for i in range(len(temp)):
  marker_genes[temp["celltype"][i]] = temp["markergenes"][i].split(",")

#calculate overlap between data derived marker genes and known marker genes
cell_annotation = sc.tl.marker_gene_overlap(adata, marker_genes, key='rank_genes_specter')
cell_annotation_norm = sc.tl.marker_gene_overlap(adata, marker_genes, key='rank_genes_specter', normalize='reference')

#annotation based on overlap
liste = []
for i in range(1,cell_annotation_norm.shape[1]+1):
	liste.append(cell_annotation_norm[str(i)].nlargest(2))

annotation_list = []
for i in range(len(liste)):
  if liste[i][0] <= 0.1: #case: no significant overlap between marker genes
    annotation_list.append("unannotated"+str(i))
  elif liste[i][0] != 0 and liste[i][1] != 0 and liste[i][1] + 0.25 >= liste[i][0] and {snakemake.params.subclustering}.pop() == "True": #case: 2 overlaps are almost equal & subclustering is active
    annotation_list.append(liste[i].index[0] + "_" + liste[i].index[1])
  elif liste[i].index[0] in annotation_list: #case: 1 or more cluster(s) are already annotated with this celltype
    count = [elem for elem in annotation_list if elem.startswith(liste[i].index[0])]
    annotation_list.append(liste[i].index[0]+str(len(count)))
  else: #case: annotate the cluster based on the highest overlap
    annotation_list.append(liste[i].index[0])

print(annotation_list)

#plot overlap as heatmap
plt.figure(figsize=(13.5,10))
plt.margins(y=0.7)
heatmap = sb.heatmap(cell_annotation_norm, cbar=False, annot=True, annot_kws={"size":15})
heatmap = heatmap.get_figure()
heatmap.savefig("figures/markergenes_overlap.png")

#define a nice colour map for gene expression
colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
colors3 = plt.cm.Greys_r(np.linspace(0.7,0.8,20))
colorsComb = np.vstack([colors3, colors2])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)

#visualizing given genes across clusters
genes = {snakemake.params.genes}.pop().split(",")

for i in range(len(genes)):
	sc.pl.umap(adata, color=genes[i], use_raw=False, color_map=mymap, save="_"+genes[i])
	sc.pl.violin(adata, groupby='specter', keys=genes[i], use_raw=False, save="_"+genes[i])

#renaming cluster categories (based on annotations)
adata.rename_categories('specter', annotation_list)

#plot annotated clustering
sc.pl.umap(adata, color='specter', size=15, legend_loc='on data', save="_annotated_clustering")

adata.write("file_dir/post_clustering_specter.h5ad")
np.save("file_dir/annotation_list_specter.npy", annotation_list)
