import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb

#read in units-file
units = pd.read_table(snakemake.params.units)

#extract columns
samples = np.asarray(units["sample"])
print(samples)
sample_strings = np.asarray(units["sample_alias"])
print(sample_strings)
region = np.asarray(units["region"])
print(region)
#check if contrast- & variable-column exist
diff_testing = False
user_col = False

if len(units.columns) == 4 or len(units.columns) == 5:
  contrast = np.asarray(units["contrast"])
  diff_testing = True
  if len(units.columns) == 5 and units[units.columns[4]].isnull().values.any() == False:
    name = units.columns[4]
    variable = np.asarray(units[name])
    user_col = True

#load data (in anndata object)
adata = sc.read("data/"+samples[0]+"/outs/filtered_feature_bc_matrix/matrix.mtx.gz") #read in matrix
adata = adata.transpose()
adata.X = adata.X.toarray()

barcodes = pd.read_csv("data/"+samples[0]+"/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", header=None, sep='\t') #read in barcodes
genes = pd.read_csv("data/"+samples[0]+"/outs/filtered_feature_bc_matrix/features.tsv.gz", header=None, sep='\t') #read in genes

#annotate data (barcodes)
barcodes.rename(columns={0:'barcode'}, inplace=True)
barcodes.set_index('barcode', inplace=True)
adata.obs = barcodes #anndata.obs contains barcodes
adata.obs['sample'] = [sample_strings[0]]*adata.n_obs
adata.obs['region'] = [region[0]]*adata.n_obs
if user_col:
  adata.obs[name] = [variable[0]]*adata.n_obs
if diff_testing:
  adata.obs["contrast"] = [contrast[0]]*adata.n_obs

#annotate data (genes)
genes.rename(columns={0:'gene_id', 1:'gene_symbol'}, inplace=True)
genes.set_index('gene_symbol', inplace=True)
adata.var = genes #anndata.var containes genes
adata.var_names_make_unique()

#adding other samples to anndata object in a for-loop
for i in range(1,len(sample_strings)):

    #Load data
    adata_tmp = sc.read("data/"+samples[i]+"/outs/filtered_feature_bc_matrix/matrix.mtx.gz")
    adata_tmp = adata_tmp.transpose()
    adata_tmp.X = adata_tmp.X.toarray()

    barcodes_tmp = pd.read_csv("data/"+samples[i]+"/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", header=None, sep='\t')
    genes_tmp = pd.read_csv("data/"+samples[i]+"/outs/filtered_feature_bc_matrix/features.tsv.gz", header=None, sep='\t')

    #Annotate data
    barcodes_tmp.rename(columns={0:'barcode'}, inplace=True)
    barcodes_tmp.set_index('barcode', inplace=True)
    adata_tmp.obs = barcodes_tmp
    adata_tmp.obs['sample'] = [sample_strings[i]]*adata_tmp.n_obs
    adata_tmp.obs['region'] = [region[i]]*adata_tmp.n_obs
    if user_col:
      adata_tmp.obs[name] = [variable[i]]*adata_tmp.n_obs
    if diff_testing:
      adata_tmp.obs["contrast"] = [contrast[i]]*adata_tmp.n_obs

    genes_tmp.rename(columns={0:'gene_id', 1:'gene_symbol'}, inplace=True)
    genes_tmp.set_index('gene_symbol', inplace=True)
    adata_tmp.var = genes_tmp
    adata_tmp.var_names_make_unique()

    adata.var_names_make_unique()
    adata = adata.concatenate(adata_tmp, batch_key='sample_id')
    adata.obs.drop(columns=['sample_id'], inplace=True)
    adata.obs_names = [c.split("-")[0] for c in adata.obs_names]
    adata.obs_names_make_unique(join='_')


print(adata.obs['region'].value_counts()) #print amount of cells per region
print('')
if user_col:
  print(adata.obs[name].value_counts()) #print amount of cells per classificator in variable-column
print('')
print(adata.obs['sample'].value_counts()) #print amount of cells per sample

adata.var.drop(columns=[2], inplace=True)

adata.write("file_dir/data.h5ad")
