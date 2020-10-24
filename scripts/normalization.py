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

#clustering for calculation of sizefactors
adata_pp = adata.copy()
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp)
sc.pp.pca(adata_pp, n_comps=15)
sc.pp.neighbors(adata_pp)
sc.tl.louvain(adata_pp, key_added='groups', resolution=0.5)

#export groups + datamatrix for sizefactors calculation in R
input_groups = adata_pp.obs['groups']
data_mat = adata.X.T
np.savetxt("file_dir/data_mat.csv", data_mat, delimiter=",")
input_groups.to_csv("file_dir/input_groups.csv")
