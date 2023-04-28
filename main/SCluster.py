###################################################################
##SClustering##
import rpy2.robjects as robjects
from rpy2.robjects import r, pandas2ri,numpy2ri
import pandas as pd
import numpy as np
import anndata
import scanpy as sc
import time

numpy2ri.activate()
pandas2ri.activate()

##导入R包
robjects.r.source('methods.R')

def SClustering(adata,c_nums=None,sc3=False,soup=False,seurat=False,RaceID3=False,
                scanpy=False,cidr=False,simlr=False,sincera=False,sharp=False,results=False):
    SC=Methods(adata,c_nums,results)
    SC.read_data()
    if sc3:
        print("Start SC3 clustering.....\n",flush=True)
        start_time=time.time()
        SC.SC_sc3()
        end_time=time.time()
        print('SC3 takes {} secs'.format(end_time-start_time),flush=True)

    if soup:
        print("Start SOUP clustering.....\n",flush=True)
        start_time=time.time()
        SC.SC_soup()
        end_time=time.time()
        print('SOUP takes {} secs'.format(end_time-start_time),flush=True)

    if cidr:
        print("Start CIDR clustering.....\n",flush=True)
        start_time=time.time()
        SC.SC_cidr()
        end_time=time.time()
        print('CIDR takes {} secs'.format(end_time-start_time),flush=True)

    if simlr:
        print("Start SIMLR clustering.....\n",flush=True)
        start_time=time.time()
        SC.SC_simlr()
        end_time=time.time()
        print('CIDR takes {} secs'.format(end_time-start_time),flush=True)

    if sincera:
        print("Start SINCERA clustering.....\n",flush=True)
        start_time=time.time()
        SC.SC_sincera()
        end_time=time.time()
        print('SINCERA takes {} secs'.format(end_time-start_time),flush=True)
    
    if sharp:
        print("Start SHARP clustering.....\n",flush=True)
        start_time=time.time()
        SC.SC_sharp()
        end_time=time.time()
        print('SHARP takes {} secs'.format(end_time-start_time),flush=True)
    
    if seurat:
        print("Start Seurat clustering.....\n",flush=True)
        start_time=time.time()
        SC.SC_seurat()
        end_time=time.time()
        print('Seurat takes {} secs'.format(end_time-start_time),flush=True)

    if RaceID3:
        print("Start RaceID3 clustering.....\n",flush=True)
        start_time=time.time()
        SC.SC_RaceID3()
        end_time=time.time()
        print('RaceID3 takes {} secs'.format(end_time-start_time),flush=True)
    
    if scanpy:
        print("Start Scanpy clustering.....\n",flush=True)
        start_time=time.time()
        SC.SC_scanpy()
        end_time=time.time()
        print('Scanpy  takes {} secs'.format(end_time-start_time),flush=True)
    
    
    return SC.adata

class Methods(object):
    def __init__(self,adata,c_nums,results):

        self.sc3_SCluster=None
        self.soup_SCluster=None
        self.simlr_SCluster=None
        self.cidr_SCluster=None
        self.sincera_SCluster=None
        self.sharp_SCluster=None
        self.seurat_SCluster=None
        self.RaceID3_SCluster=None

        self.results=results
        self.sc3out=None
        self.soupout=None
        self.simlrout=None
        self.cidrout=None
        self.sinceraout=None
        self.sharpout=None
        self.seuratout=None
        self.RaceID3out=None

        self.matrix=None
        self.c_nums=c_nums
        self.adata=adata
        self.results=results

    def SC_sc3(self):
        robjects.r.library("SC3")
        robjects.r.library("SingleCellExperiment")
        self.sc3_SCluster=robjects.r["sc3_SCluster"]
        self.sc3out=self.sc3_SCluster(self.c_nums,self.matrix,self.results)
        self.adata.obs['sc3']=self.sc3out[-1]
        self.adata.obs['sc3']=self.adata.obs['sc3'].astype('category')

    def SC_soup(self):
        robjects.r.library("SOUP")
        self.soup_SCluster=robjects.r["soup_SCluster"]
        self.soupout=self.soup_SCluster(self.c_nums,self.matrix,self.results)
        self.adata.obs['soup']=self.soupout[-1]
        self.adata.obs['soup']=self.adata.obs['soup'].astype('category')
    
    def SC_simlr(self):
        robjects.r.library("SIMLR")
        self.simlr_SCluster=robjects.r["simlr_SCluster"]
        self.simlrout=self.simlr_SCluster(self.c_nums,self.matrix,self.results)
        self.adata.obs['simlr']=self.simlrout[-1]
        self.adata.obs['simlr']=self.adata.obs['simlr'].astype('category')

    
    def SC_cidr(self):
        robjects.r.library("cidr")
        self.cidr_SCluster=robjects.r["cidr_SCluster"]
        self.cidrout=self.cidr_SCluster(self.c_nums,self.matrix,self.results)
        self.adata.obs['cidr']=self.cidrout[-1]
        self.adata.obs['cidr']=self.adata.obs['cidr'].astype('category')
    
    def SC_sincera(self):
        robjects.r.library("SINCERA")
        robjects.r.library("cluster")
        self.sincera_SCluster=robjects.r["sincera_SCluster"]
        self.sinceraout=self.sincera_SCluster(self.c_nums,self.matrix,self.results)
        self.adata.obs['sincera']=self.sinceraout[-1]
        self.adata.obs['sincera']=self.adata.obs['sincera'].astype('category')

    def SC_sharp(self):
        robjects.r.library("SHARP")
        self.sharp_SCluster=robjects.r["sharp_SCluster"]
        self.sharpout=self.sharp_SCluster(self.c_nums,self.matrix,self.results)
        self.adata.obs['sharp']=self.sharpout[-1]
        self.adata.obs['sharp']=self.adata.obs['sharp'].astype('category')

    def SC_seurat(self):
        robjects.r.library("Seurat")
        self.seurat_SCluster=robjects.r["seurat_SCluster"]
        self.seuratout=self.seurat_SCluster(self.c_nums,self.matrix,self.results)
        self.adata.obs['seurat']=self.seuratout[-1]
        self.adata.obs['seurat']=self.adata.obs['seurat'].astype('category')
    
    def SC_RaceID3(self):
        robjects.r.library("RaceID")
        self.RaceID3_SCluster=robjects.r["RaceID3_SCluster"]
        self.RaceID3out=self.RaceID3_SCluster(self.c_nums,self.matrix,self.results)
        self.adata.obs['RaceID3']=self.RaceID3out[-1]
        self.adata.obs['RaceID3']=self.adata.obs['RaceID3'].astype('category')
    
    def SC_scanpy(self):
        self.scanpyout=self.adata.copy()
        self.scanpyout.var_names_make_unique()
        sc.pp.normalize_total(self.scanpyout, target_sum=1e4)
        sc.pp.log1p(self.scanpyout)
        sc.pp.highly_variable_genes(self.scanpyout, min_mean=0.0125, max_mean=3, min_disp=0.5)
        self.scanpyout = self.scanpyout[:, self.scanpyout.var.highly_variable]
        sc.tl.pca(self.scanpyout, svd_solver='arpack')
        sc.pp.neighbors(self.scanpyout)
        for i in range(1,101):
            sc.tl.leiden(self.scanpyout,resolution=(i)*0.02)
            cluster=pd.DataFrame(self.scanpyout.obs["leiden"])['leiden']
            if len(np.unique(cluster))==int(self.c_nums):
                break
            if i==30:
                sc.tl.leiden(self.scanpyout)
        if self.results:
            self.scanpyout.write("scanpy_SCluster.h5ad",compression='gzip')
        self.adata.obs['scanpy']=self.scanpyout.obs['leiden']
        self.adata.obs['scanpy']=self.adata.obs['scanpy'].astype('category')
    
    def read_data(self):
        self.matrix=pd.DataFrame(self.adata.X.todense())
        self.matrix.index=self.adata.obs.index
        self.matrix.columns=self.adata.var.index
        self.matrix=self.matrix.T
        print('{} samples {} features'.format(self.matrix.shape[1],self.matrix.shape[0]),flush=True)
    










