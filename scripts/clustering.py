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

#perform clustering - using highly variable genes
sc.tl.louvain(adata, key_added='louvain_r1')
sc.tl.louvain(adata, resolution=float({snakemake.params.clustering_resolution}.pop()), key_added='louvain_r'+str(snakemake.params.clustering_resolution), random_state=10)

#visualize the clustering and how this is reflected by different technical covariates
sc.pl.umap(adata, color=['louvain_r1', 'louvain_r'+str(snakemake.params.clustering_resolution)], palette=sc.pl.palettes.vega_20, save="_clustering_result")
sc.pl.umap(adata, color=['region', 'n_counts'], save="_region_ncounts_oncluster")
sc.pl.umap(adata, color=['log_counts', 'mt_frac'], save="_logcounts_mtfrac_oncluster")

#marker genes & cluster annotation
#calculate marker genes derived from data
sc.tl.rank_genes_groups(adata, groupby='louvain_r'+str(snakemake.params.clustering_resolution), key_added='rank_genes_r'+str(snakemake.params.clustering_resolution))

#plot marker genes
sc.pl.rank_genes_groups(adata, key='rank_genes_r'+str(snakemake.params.clustering_resolution), fontsize=12, save="_markergenes")

#read in specified  marker genes
temp = pd.read_table({snakemake.params.celltype_annotation}.pop())
marker_genes = {}
for i in range(len(temp)):
  marker_genes[temp["celltype"][i]] = temp["markergenes"][i].split(",")

#calculate overlap between data derived marker genes and known marker genes
cell_annotation = sc.tl.marker_gene_overlap(adata, marker_genes, key='rank_genes_r'+str(snakemake.params.clustering_resolution))
cell_annotation_norm = sc.tl.marker_gene_overlap(adata, marker_genes, key='rank_genes_r'+str(snakemake.params.clustering_resolution), normalize='reference')

#annotation based on overlap
liste = []
for i in range(cell_annotation_norm.shape[1]):
	liste.append(cell_annotation_norm[str(i)].nlargest(2))

annotation_list = []
for i in range(len(liste)):
  if liste[i][0] <= 0.1: #case: no significant overlap between marker genes
    annotation_list.append("unannotated"+str(i))
  elif liste[i][0] != 0 and liste[i][1] != 0 and liste[i][1] + 0.25 >= liste[i][0] and {snakemake.params.subclustering}.pop() == "True": #case: 2 overlaps are almost equal & subclustering is activ in config-file
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

#visualizing specified genes across clusters
genes = {snakemake.params.genes}.pop().split(",")

for i in range(len(genes)):
	sc.pl.umap(adata, color=genes[i], use_raw=False, color_map=mymap, save="_"+genes[i])
	sc.pl.violin(adata, groupby='louvain_r'+str(snakemake.params.clustering_resolution), keys=genes[i], use_raw=False, save="_"+genes[i])

#renaming cluster categories (based on annotations)
adata.rename_categories('louvain_r'+str(snakemake.params.clustering_resolution), annotation_list)

#plot annotated clustering
sc.pl.umap(adata, color='louvain_r'+str(snakemake.params.clustering_resolution), size=15, legend_loc='on data', save="_annotated_clustering")

adata.write("file_dir/post_clustering_louvain.h5ad")
np.save("file_dir/annotation_list_louvain.npy", annotation_list)
