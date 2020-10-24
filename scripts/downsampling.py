import scanpy as sc
import numpy as np
#import scipy as sp
import pandas as pd

adata = sc.read_h5ad(snakemake.input[0])
down_sample = np.loadtxt(snakemake.input[1], delimiter=",")

adata.obs["subset"] = down_sample
print(adata.obs)
print(adata)

adata = adata[adata.obs["subset"] == 1]
print(adata)
print(adata.obs)

adata.obs.drop(labels="subset", axis=1, inplace=True)
print(adata.obs)

adata.write("file_dir/normed_data_sphetcher.h5ad")
