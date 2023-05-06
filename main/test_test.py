from  SCluster import SCluster
import scanpy as sc
import anndata
import pandas as pd
import anndata as ad


data=pd.read_csv('/nfs_genome1/wanxinjiang/single_cell_benchmark/data/Metadata_Esophagus/Matrix_Esophagus.txt.gz',sep='\t')
type=pd.read_csv('/nfs_genome1/wanxinjiang/single_cell_benchmark/data/Metadata_Esophagus/celltype.csv')
type=type['type']
adata=ad.AnnData(data.T)
adata.obs["annotation"]=type.values


sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

adata.raw=adata

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata= adata[:, adata.var.highly_variable]

SCluster_res,adata_res=SCluster(adata=adata,c_nums=11,\
                           sc3=True,\
                           simlr=True,scanpy=True,\
                           soup=True,sharp=True,\
                           first_result=True)

SCluster_res.write('/nfs_genome1/wanxinjiang/SCluster/demo/Esophagus_SCluster.h5ad',compression='gzip')
adata_res.write('/nfs_genome1/wanxinjiang/SCluster/demo/Esophagus_res.h5ad',compression='gzip')
