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

#read in sizefactors and add them to anndata.obs
size_factors = pd.read_csv(snakemake.input[1])
size_factors_array = size_factors.to_numpy()
adata.obs["size_factors"] = size_factors_array

#plotting
sc.pl.scatter(adata, 'size_factors', 'n_counts', save="_normed_ncounts_sizefactors", show=False) #Plotten
sc.pl.scatter(adata, 'size_factors', 'n_genes', save="_normed_ngenes_sizefactors", show=False)

plt.figure(1)
distplot = sb.distplot(size_factors, bins=50, kde=False)
distplot = distplot.get_figure()
distplot.savefig("figures/sizefactors_plot.png")

#save the count data in a counts layer
adata.layers["counts"] = adata.X.copy()
print(adata.layers["counts"])

#normalize adata
adata.X /= adata.obs['size_factors'].values[:,None]
sc.pp.log1p(adata)

#store the full data set in 'raw' as log-normalised data for statistical testing
adata.X = sp.sparse.csr_matrix(adata.X)
adata.raw = adata

#batch_Correction
sc.pp.combat(adata, key="sample")

#compute highly variable genes(HVGs)
sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=int(snakemake.params.amount_of_hvgs))
print('\n','Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))

#plot (HVGs)
sc.pl.highly_variable_genes(adata, save="_high_variable_genes")

#prepare anndata object for downsampling with sphetcher if specified in config-file
if snakemake.params.downsampling == "sphetcher":
  sc.pp.pca(adata, n_comps=50, use_highly_variable=False, svd_solver='arpack')
  sc.pp.neighbors(adata)
  np.savetxt("file_dir/sphetcher_pca.csv", adata.obsm["X_pca"], delimiter=",")

adata.write("file_dir/normed_data_normal.h5ad")
