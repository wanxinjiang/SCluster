
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
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import fowlkes_mallows_score
```
### Import data


