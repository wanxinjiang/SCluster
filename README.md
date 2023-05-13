
# SClustering
SCluster: An integrated clustering framework for single-cell RNA-seq data
## Installation
```
mkdir single_cell
tar -xzvf single_cell.tar.gz  -C ./single_cell/
```
## SClustering Examples
### Import python modules
```
import anndata
import anndata as ad
import pandas as pd 
import numpy as np
import seaborn as sns
import scanpy as sc
from SCluster import SCluster
import matplotlib.pyplot as plt
from sklearn.metrics.cluster import adjusted_rand_score
```
### Import data
```
adata=anndata.read('./datasets/Trachea.h5ad')
>>>adata
AnnData object with n_obs × n_vars = 4798 × 17197
    obs: 'annotation'
```
### Filtering cells and genes
```
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
```
```
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
```
Actually do the filtering by slicing the AnnData object.
```
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
```
```
adata.raw=adata
```
Total-count normalize (library-size correct) the data matrix X to 10,000 reads per cell, so that counts become comparable among cells.
Logarithmize the data.
Identify highly-variable genes.
Actually do the filtering.
```
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata= adata[:, adata.var.highly_variable]
>>>adata
View of AnnData object with n_obs × n_vars = 3819 × 3256
    obs: 'annotation', 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt'
    var: 'n_cells', 'mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
    uns: 'log1p', 'hvg'
```
### Performing SCluster clustering
```
adata=SCluster(adata=adata,\
               sc3=True,cidr=True,sharp=True,scanpy=True,\
               soup=True,seurat=True,simlr=True,RaceID3=True,sincera=True)
>>>adata
AnnData object with n_obs × n_vars = 3819 × 3256
    obs: 'annotation', 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'sc3', 'soup', 'cidr', 'sincera', 'sharp', 'scanpy', 'seurat', 'RaceID3', 'simlr', 'SCluster'
    var: 'n_cells', 'mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
    uns: 'log1p', 'hvg', 'neighbors', 'umap'
    obsm: 'X_pca', 'X_umap'
    obsp: 'distances', 'connectivities'
```
```
>>>adjusted_rand_score(adata.obs['annotation'],adata.obs['SCluster'])
0.9688839138590054
```
### Plot the clusters
```
plt.rcParams['figure.dpi']=300
sc.set_figure_params(fontsize=10)
adata.obs['annotation']=adata.obs['annotation'].astype('category')
sc.pl.umap(adata,color=['annotation','SCluster','soup','sharp','scanpy','seurat','sc3','sincera','cidr','simlr','RaceID3'],legend_fontsize=8,s=12,frameon=False)
```
![这是图片](/Examples/umap.png "umap")
