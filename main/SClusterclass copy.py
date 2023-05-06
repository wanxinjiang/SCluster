


###################################################################
##SClustering##
import rpy2.robjects as robjects
from rpy2.robjects import r, pandas2ri,numpy2ri
import pandas as pd
import numpy as np
import anndata
import scanpy as sc
import time
import anndata as ad
from multiprocessing import Pool,Lock,Process

numpy2ri.activate()
pandas2ri.activate()

# robjects.r('options(warn=-1)')
# robjects.r('options(message=-1)')
# lock=Lock()


##导入R包
robjects.r.source('methods.R')

def SClustering(adata,c_nums=None,sc3=False,soup=False,seurat=False,RaceID3=False,
                scanpy=False,cidr=False,simlr=False,sincera=False,sharp=False,results=False):
    
    SC=Methods(adata,c_nums,results)
    SC.read_data()
    rc=run_clustering(SC)
    
    
    pool = Pool(20)

    pool.apply_async(rc.run_sc3,(sc3,))
    pool.apply_async(rc.run_soup,(soup,))
    pool.apply_async(rc.run_cidr,(cidr,))
    pool.apply_async(rc.run_simlr,(simlr,))
    pool.apply_async(rc.run_sincera,(sincera,))
    pool.apply_async(rc.run_sharp,(sharp,))
    pool.apply_async(rc.run_RaceID3,(RaceID3,))
    pool.apply_async(rc.run_seurat,(seurat,))
    pool.apply_async(rc.run_scanpy,(scanpy,))
    
    pool.close()
    pool.join()
    
    run_results=rc.run_results
    print(run_results)
    for result in run_results:
        if result!=None:
            adata.obs[result[0]]=result[1].values

    
    return adata

class run_clustering(object):
    def __init__(self,SC):
        self.run_results=[]
        self.SC=SC
    
    def run_sc3(self,sc3):
        if sc3:
            print("Start SC3 clustering...",flush=True)
            self.SC.SC_sc3()
            print('SC3 done.')
            self.run_results.append(['sc3',self.SC.adata.obs["sc3"]])

    def run_cidr(self,cidr):
        if cidr:
            print("Start CIDR clustering...",flush=True)
            self.SC.SC_cidr()
            print('CIDR done.')
            self.run_results.append(['cidr',self.SC.adata.obs["cidr"]])

    def run_soup(self,soup):
        if soup:
            print("Start SOUP clustering...",flush=True)
            self.SC.SC_soup()
            print('SOUP done.')
            self.run_results.append(['soup',self.SC.adata.obs["soup"]])

    def run_simlr(self,simlr):
        if simlr:
            print("Start SIMLR clustering...",flush=True)
            self.SC.SC_simlr()
            print('SIMLR done.')
            self.run_results.append(['simlr',self.SC.adata.obs["simlr"]])
            
    def run_RaceID3(self,RaceID3):
        if RaceID3:
            print("Start RaceID3 clustering...",flush=True)
            self.SC.SC_RaceID3()
            print('RaceID3 done.')
            self.run_results.append(['RaceID3',self.SC.adata.obs["RaceID3"]])

    def run_seurat(self,seurat):
        if seurat:
            print("Start Seurat clustering...",flush=True)
            self.SC.SC_seurat()
            print('Seurat done.')
            self.run_results.append(['seurat',self.SC.adata.obs["seurat"]])

    def run_sincera(self,sincera):
        if sincera:
            print("Start SINCERA clustering...",flush=True)
            self.SC.SC_sincera()
            print('SINCERA done.')
            self.run_results.append(['sincera',self.SC.adata.obs["sincera"]])

    def run_sharp(self,sharp):
        if sharp:
            print("Start SHARP clustering...",flush=True)
            self.SC.SC_sharp()
            print('SHARP done.')
            self.run_results.append(['sharp',self.SC.adata.obs["sharp"]])

    def run_scanpy(self,scanpy):
        if scanpy:
            print("Start Scanpy  clustering...",flush=True)
            self.SC.SC_scanpy()
            print('Scanpy done.')
            self.run_results.append(['scanpy',self.SC.adata.obs["scanpy"]])
        

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
        self.matrix_raw=None
        self.c_nums=c_nums
        self.adata=adata
        self.results=results
        self.adata_raw=adata.raw

    def SC_sc3(self):
        robjects.r.library("SC3")
        robjects.r.library("SingleCellExperiment")
        self.sc3_SCluster=robjects.r["sc3_SCluster"]
        self.sc3out=self.sc3_SCluster(self.c_nums,self.matrix,self.results)
        # lock.acquire()
        self.adata.obs['sc3']=self.sc3out[-1]
        self.adata.obs['sc3']=self.adata.obs['sc3'].astype('category')
        # lock.release()

    def SC_soup(self):
        robjects.r.library("SOUP")
        self.soup_SCluster=robjects.r["soup_SCluster"]
        self.soupout=self.soup_SCluster(self.c_nums,self.matrix,self.results)
        # lock.acquire()
        self.adata.obs['soup']=self.soupout[-1]
        self.adata.obs['soup']=self.adata.obs['soup'].astype('category')
        # lock.release()
    
    def SC_simlr(self):
        robjects.r.library("SIMLR")
        self.simlr_SCluster=robjects.r["simlr_SCluster"]
        self.simlrout=self.simlr_SCluster(self.c_nums,self.matrix,self.results)
        # lock.acquire()
        self.adata.obs['simlr']=self.simlrout[-1]
        self.adata.obs['simlr']=self.adata.obs['simlr'].astype('category')
        # lock.release()

    
    def SC_cidr(self):
        robjects.r.library("cidr")
        self.cidr_SCluster=robjects.r["cidr_SCluster"]
        self.cidrout=self.cidr_SCluster(self.c_nums,self.matrix_raw,self.results)
        # lock.acquire()
        self.adata.obs['cidr']=self.cidrout[-1]
        self.adata.obs['cidr']=self.adata.obs['cidr'].astype('category')
        # lock.release()
    
    def SC_sincera(self):
        robjects.r.library("SINCERA")
        robjects.r.library("cluster")
        self.sincera_SCluster=robjects.r["sincera_SCluster"]
        self.sinceraout=self.sincera_SCluster(self.c_nums,self.matrix_raw,self.results)
        # lock.acquire()
        self.adata.obs['sincera']=self.sinceraout[-1]
        self.adata.obs['sincera']=self.adata.obs['sincera'].astype('category')
        # lock.release()

    def SC_sharp(self):
        robjects.r.library("SHARP")
        self.sharp_SCluster=robjects.r["sharp_SCluster"]
        self.sharpout=self.sharp_SCluster(self.c_nums,self.matrix,self.results)
        # lock.acquire()
        self.adata.obs['sharp']=self.sharpout[-1]
        self.adata.obs['sharp']=self.adata.obs['sharp'].astype('category')
        # lock.release()

    def SC_seurat(self):
        robjects.r.library("Seurat")
        self.seurat_SCluster=robjects.r["seurat_SCluster"]
        self.seuratout=self.seurat_SCluster(self.c_nums,self.matrix,self.results)
        # lock.acquire()
        self.adata.obs['seurat']=self.seuratout[-1]
        self.adata.obs['seurat']=self.adata.obs['seurat'].astype('category')
        # lock.release()
    

    def SC_RaceID3(self):
        robjects.r.library("RaceID")
        self.RaceID3_SCluster=robjects.r["RaceID3_SCluster"]
        self.RaceID3out=self.RaceID3_SCluster(self.c_nums,self.matrix_raw,self.results)
        # lock.acquire()
        self.adata.obs['RaceID3']=self.RaceID3out[-1]
        self.adata.obs['RaceID3']=self.adata.obs['RaceID3'].astype('category')
        # lock.release()
    
    def SC_scanpy(self):
        self.scanpyout=self.adata.copy()
        self.scanpyout.var_names_make_unique()
        sc.tl.pca(self.scanpyout, svd_solver='arpack')
        sc.pp.neighbors(self.scanpyout)
        sc.tl.leiden(self.scanpyout)
        if self.results:
            self.scanpyout.write("scanpy_SCluster.h5ad",compression='gzip')
        # lock.acquire()
        self.adata.obs['scanpy']=self.scanpyout.obs['leiden']
        self.adata.obs['scanpy']=self.adata.obs['scanpy'].astype('category')
        # lock.release()
    
    def read_data(self):
        try:
            self.matrix=pd.DataFrame(self.adata.X.todense())
        except:
            self.matrix=pd.DataFrame(self.adata.X)
            
        try:
            self.matrix_raw=pd.DataFrame(self.adata_raw.X.todense())
        except:
            self.matrix_raw=pd.DataFrame(self.adata_raw.X)
            
        self.matrix.index=self.adata.obs.index
        self.matrix.columns=self.adata.var.index
        self.matrix=self.matrix.T
        
        self.matrix_raw.index=self.adata.obs.index
        self.matrix_raw.columns=self.adata_raw.var.index
        self.matrix_raw=self.matrix_raw.T
        print('Raw expression',flush=True)
        print('{} samples {} features\n'.format(self.matrix_raw.shape[1],self.matrix_raw.shape[0]),flush=True)
        print('After selecting highly variable',flush=True)
        print('{} samples {} features\n'.format(self.matrix.shape[1],self.matrix.shape[0]),flush=True)
    

    







