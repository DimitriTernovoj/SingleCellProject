import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb

adata = sc.read_h5ad(snakemake.input[0])

#compute variables that are interesting for quality control
adata.obs['n_counts'] = adata.X.sum(1)
adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
adata.obs['n_genes'] = (adata.X > 0).sum(1)

mt_gene_mask = [gene.startswith('mt-') for gene in adata.var_names]
adata.obs['mt_frac'] = adata.X[:, mt_gene_mask].sum(1)/adata.obs['n_counts']

#visualizations
sc.pl.violin(adata, 'n_counts', groupby='sample', size=2, log=True, cut=0, save="_ncounts_groupedby_sample", show=False)
sc.pl.violin(adata, 'mt_frac', groupby='sample', save="_mt_frac_groupedby_sample", show=False)
sc.pl.scatter(adata, 'n_counts', 'n_genes', color='mt_frac', save="_ncounts_ngenes", show=False)

#specifieng cut-offs
counts = adata.obs["n_counts"].to_numpy()
gene_counts = adata.obs["n_genes"].to_numpy()

if snakemake.params.lower_quantile_counts != "":
  min_counts_cells = round(np.quantile(counts,float(snakemake.params.lower_quantile_counts)))
else:
  min_counts_cells = int(snakemake.params.min_counts)

if snakemake.params.upper_quantile_counts != "":
  max_counts_cells = round(np.quantile(counts,float(snakemake.params.upper_quantile_counts)))
else:
  max_counts_cells = int(snakemake.params.max_counts)

if snakemake.params.lower_quantile_genes != "":
  min_counts_genes = round(np.quantile(gene_counts,float(snakemake.params.lower_quantile_genes)))
else:
  min_counts_genes = int(snakemake.params.min_genes)

#further visualizations
plt.figure(1)
p3 = sb.distplot(adata.obs['n_counts'], kde=False)
plt.axvline(x = min_counts_cells, color="red")
plt.axvline(x = max_counts_cells, color="red")
p3_fig = p3.get_figure()
p3_fig.savefig("figures/n_counts.png")

plt.figure(2)
p4 = sb.distplot(adata.obs['n_counts'][adata.obs['n_counts']<np.quantile(adata.obs["n_counts"], 0.5)], kde=False, bins=60)
plt.axvline(x = min_counts_cells, color="red")
p4_fig = p4.get_figure()
p4_fig.savefig("figures/lower_ncounts.png")

plt.figure(3)
p5 = sb.distplot(adata.obs['n_counts'][adata.obs['n_counts']>np.quantile(adata.obs["n_counts"], 0.5)], kde=False, bins=60)
plt.axvline(x = max_counts_cells, color="red")
p5_fig = p5.get_figure()
p5_fig.savefig("figures/upper_ncounts.png")

plt.figure(4)
p6 = sb.distplot(adata.obs['n_genes'], kde=False, bins=60)
plt.axvline(x = min_counts_genes, color="red")
p6_fig = p6.get_figure()
p6_fig.savefig("figures/ncounts_genes.png")

plt.figure(5)
p7 = sb.distplot(adata.obs['n_genes'][adata.obs['n_genes']<np.quantile(adata.obs["n_genes"], 0.5)], kde=False, bins=60)
plt.axvline(x = min_counts_genes, color="red")
p7_fig = p7.get_figure()
p7_fig.savefig("figures/lower_ncounts_genes.png")


#print information after qc on cells
print('Total number of cells: {:d}'.format(adata.n_obs))

sc.pp.filter_cells(adata, min_counts = min_counts_cells)
print('Number of cells after min count filter: {:d}'.format(adata.n_obs))

sc.pp.filter_cells(adata, max_counts = max_counts_cells)
print('Number of cells after max count filter: {:d}'.format(adata.n_obs))

adata = adata[adata.obs['mt_frac'] < float(snakemake.params.mt_frac)]
print('Number of cells after MT filter: {:d}'.format(adata.n_obs))

sc.pp.filter_cells(adata, min_genes = min_counts_genes)
print('Number of cells after gene filter: {:d}'.format(adata.n_obs))

#filter genes & print number of cells after qc
sc.pp.filter_genes(adata, min_cells=int(snakemake.params.min_cells))
print('Number of genes after cell filter: {:d}'.format(adata.n_vars))

adata.write("file_dir/filtered_data.h5ad")
