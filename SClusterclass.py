##SClustering##
import rpy2.robjects as robjects
from rpy2.robjects import r, pandas2ri,numpy2ri
import pandas as pd
import numpy as np
import anndata
import scanpy as sc
import time
import os
import math
import subprocess
import anndata as ad
from multiprocessing import Pool,Lock,Process
from sklearn.metrics.cluster import davies_bouldin_score,calinski_harabasz_score
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger

import numba
numba.config.THREADING_LAYER = 'workqueue'

numpy2ri.activate()
pandas2ri.activate()


rpy2_logger.setLevel(50)

_temp='./temp/'
_seuratout_csv=_temp+'seuratout.csv'
_RaceID3out_csv=_temp+'RaceID3out.csv'
_simlrout_csv=_temp+'simlrout.csv'
_matrix_raw=_temp+'matrix_raw.csv'
_matrix=_temp+'matrix.csv'


# lock=Lock()
robjects.r.source("./SCluster/methods.R")

def SClustering(adata,c_nums=None,sc3=False,soup=False,seurat=False,RaceID3=False,
                scanpy=False,cidr=False,simlr=False,sincera=False,sharp=False,results=False):
    
    if c_nums==None:
        
        print('Start estimate number of clusters ...',flush=True)
        try:
            matrix=adata.X.todense()
        except:
            matrix=adata.X
            
        scanpyout=adata.copy()
        scanpyout.var_names_make_unique()
        sc.tl.pca(scanpyout, svd_solver='arpack')
        sc.pp.neighbors(scanpyout)
        # Estimate_K_db={}
        # Estimate_K_ch={}
        Estimate_K={}
        for i in range(1,41):
                sc.tl.leiden(scanpyout,resolution=(i)*0.1)
                cluster=np.array(scanpyout.obs["leiden"])
                if len(np.unique(cluster))>1:
                    db=davies_bouldin_score(np.array(matrix),cluster)
                    ch=calinski_harabasz_score(np.array(matrix),cluster)
                    Estimate_K[len(np.unique(cluster))]=np.log10(ch)+np.log10(db)
                    # Estimate_K_db[len(np.unique(cluster))]=np.log10(db)
                    # Estimate_K_ch[len(np.unique(cluster))]=np.log10(ch)
        c_nums=int(list(Estimate_K.keys())[list(Estimate_K.values()).index(max(Estimate_K.values()))])
        # c_nums_db=int(list(Estimate_K_db.keys())[list(Estimate_K_db.values()).index(min(Estimate_K_db.values()))])
        # c_nums_ch=int(list(Estimate_K_ch.keys())[list(Estimate_K_ch.values()).index(max(Estimate_K_ch.values()))])
        # print('db :{} ch:{}'.format(c_nums_db,c_nums_ch))
        # c_nums=math.ceil((c_nums_ch+c_nums_db)/2)
        
    print('The number of Clusters is {}.'.format(c_nums),flush=True)
    
    SC=Methods(adata,c_nums,results)
    
    ##Seurat、RacID3、SIMLR## 
    if simlr:
        print("Start SIMLR clustering...",flush=True)
        simlr_args = ["Rscript",'--vanilla', "/nfs_genome1/wanxinjiang/SCluster/SCluster/Rscript/run_simlr.R",]+[_simlrout_csv,str(c_nums),_matrix]
        subprocess.Popen(simlr_args)
        
    if seurat:
        print("Start Seurat clustering...",flush=True)
        seurat_args = ["Rscript",'--vanilla', "/nfs_genome1/wanxinjiang/SCluster/SCluster/Rscript/run_seurat.R",]+[_seuratout_csv,str(c_nums),_matrix_raw]
        subprocess.Popen(seurat_args)
        
    if RaceID3:
        print("Start RaceID3 clustering...",flush=True)
        RaceID3_args = ["Rscript",'--vanilla', "/nfs_genome1/wanxinjiang/SCluster/SCluster/Rscript/run_RaceID3.R",]+[_RaceID3out_csv,str(c_nums),_matrix_raw]
        subprocess.Popen(RaceID3_args)
    
    multi_res=[]
    run_results=[]
    
    pool = Pool(25)
    multi_res.append(pool.apply_async(run_sc3,(SC,sc3,)))
    multi_res.append(pool.apply_async(run_soup,(SC,soup,)))
    multi_res.append(pool.apply_async(run_cidr,(SC,cidr,)))
    multi_res.append(pool.apply_async(run_sincera,(SC,sincera,)))
    multi_res.append(pool.apply_async(run_sharp,(SC,sharp,)))
    multi_res.append(pool.apply_async(run_scanpy,(SC,scanpy,)))
    
    run_results=[res.get() for res in multi_res]
    
    print('*'*100)
    ssr=len([1 for m in [seurat,simlr,RaceID3] if m ])

    print('$'*100)
    ssr_flag=0
    
    while True:
        sr=[i[0] for i in run_results if i!=None] 
        
        if seurat and os.path.exists(_seuratout_csv) and 'seurat' not in sr:
            run_results.append(run_seurat(seurat))
            os.remove(_seuratout_csv)
            ssr_flag=ssr_flag+1
            
        if RaceID3 and os.path.exists(_RaceID3out_csv) and 'RaceID3' not in sr:
            run_results.append(run_RaceID3(RaceID3))
            os.remove(_RaceID3out_csv)
            ssr_flag=ssr_flag+1
            
        if simlr and os.path.exists(_simlrout_csv) and 'simlr' not in sr:
            run_results.append(run_simlr(simlr))
            os.remove(_simlrout_csv)
            ssr_flag=ssr_flag+1
            
        if ssr_flag ==ssr:
            os.remove(_matrix_raw)
            os.remove(_matrix)
            break
        
    for result in run_results:
        if result!=None:
            try:
                adata.obs[result[0]]=result[1].values
                adata.obs[result[0]]=adata.obs[result[0]].astype('category')
            except:
                adata.obs[result[0]]=result[1]
                adata.obs[result[0]]=adata.obs[result[0]].astype('category')
    
    return adata


           
def run_sc3(SC,sc3):
    if sc3:
        print("Start SC3 clustering...",flush=True)
        SC.SC_sc3()
        print('SC3 done.')
        return ['sc3',SC.adata.obs["sc3"]]

def run_cidr(SC,cidr):
    if cidr:
        print("Start CIDR clustering...",flush=True)
        SC.SC_cidr()
        print('CIDR done.')
        return ['cidr',SC.adata.obs["cidr"]]

def run_soup(SC,soup):
    if soup:
        print("Start SOUP clustering...",flush=True)
        SC.SC_soup()
        print('SOUP done.')
        return ['soup',SC.adata.obs["soup"]]

def run_simlr(simlr):
    if simlr:
        simlrout=pd.read_csv(_simlrout_csv)
        if simlrout.shape[0]>10:
            simlrout=simlrout.iloc[:,-1:].values
            print('SIMLR done.')
            return ['simlr',simlrout]
        else:
            print('SIMLR erro.')
            return None

def run_RaceID3(RaceID3):
    if RaceID3:
        RaceID3out=pd.read_csv(_RaceID3out_csv)
        if RaceID3out.shape[0]>10:
            RaceID3out=RaceID3out.iloc[:,-1:].values
            print('RaceID3 done.')
            return ['RaceID3',RaceID3out]
        else:
            print('RaceID3 erro.')
            return None

def run_seurat(seurat):
    if seurat:
        print('Seurat done.')
        seuratout=pd.read_csv(_seuratout_csv)
        if seuratout.shape[0]>10:
            seuratout=seuratout.iloc[:,-1:].values
            return ['seurat',seuratout]
        else:
            return None


def run_sincera(SC,sincera):
    if sincera:
        print("Start SINCERA clustering...",flush=True)
        SC.SC_sincera()
        print('SINCERA done.')
        return ['sincera',SC.adata.obs["sincera"]]

def run_sharp(SC,sharp):
    if sharp:
        print("Start SHARP clustering...",flush=True)
        SC.SC_sharp()
        print('SHARP done.')
        return ['sharp',SC.adata.obs["sharp"]]

def run_scanpy(SC,scanpy):
    if scanpy:
        print("Start Scanpy  clustering...",flush=True)
        SC.SC_scanpy()
        print('Scanpy done.')
        return ['scanpy',SC.adata.obs["scanpy"]]
        

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
        
        self.matrix.to_csv(_matrix)
        self.matrix_raw.to_csv(_matrix_raw)
        
        print('Raw expression',flush=True)
        print('{} samples {} features\n'.format(self.matrix_raw.shape[1],self.matrix_raw.shape[0]),flush=True)
        print('After selecting highly variable',flush=True)
        print('{} samples {} features\n'.format(self.matrix.shape[1],self.matrix.shape[0]),flush=True)

    
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
    
    # def SC_simlr(self):
    #     robjects.r.library("SIMLR")
    #     self.simlr_SCluster=robjects.r["simlr_SCluster"]
    #     self.simlrout=self.simlr_SCluster(self.c_nums,self.matrix,self.results)
    #     # lock.acquire()
    #     self.adata.obs['simlr']=self.simlrout[-1]
    #     self.adata.obs['simlr']=self.adata.obs['simlr'].astype('category')
    #     # lock.release()

    
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

    # def SC_seurat(self):
    #     robjects.r.library("Seurat")
    #     self.seurat_SCluster=robjects.r["seurat_SCluster"]
    #     self.seuratout=self.seurat_SCluster(self.c_nums,self.matrix_raw,self.results)
    #     # lock.acquire()
    #     self.adata.obs['seurat']=self.seuratout[-1]
    #     self.adata.obs['seurat']=self.adata.obs['seurat'].astype('category')
    #     # lock.release()
    

    # def SC_RaceID3(self):
    #     robjects.r.library("RaceID")
    #     self.RaceID3_SCluster=robjects.r["RaceID3_SCluster"]
    #     self.RaceID3out=self.RaceID3_SCluster(self.c_nums,self.matrix_raw,self.results)
    #     # lock.acquire()
    #     self.adata.obs['RaceID3']=self.RaceID3out[-1]
    #     self.adata.obs['RaceID3']=self.adata.obs['RaceID3'].astype('category')
    #     # lock.release()
    
    def SC_scanpy(self):
        self.scanpyout=self.adata.copy()
        self.scanpyout.var_names_make_unique()
        sc.tl.pca(self.scanpyout, svd_solver='arpack')
        sc.pp.neighbors(self.scanpyout)
        if int(self.c_nums)<=3:
            sc.tl.leiden(self.scanpyout,resolution=0.01)
            cluster=self.scanpyout.obs["leiden"]  
        else:
            for i in range(1,101):
                sc.tl.leiden(self.scanpyout,resolution=(i)*0.02)
                cluster=self.scanpyout.obs["leiden"]   
                if abs(len(np.unique(cluster))==int(self.c_nums)):
                            break     
                if i==101:
                    sc.tl.leiden(self.scanpyout)
                    cluster=self.scanpyout.obs["leiden"]   

        # lock.acquire()
        self.adata.obs['scanpy']=cluster
        self.adata.obs['scanpy']=self.adata.obs['scanpy'].astype('category')
        # lock.release()
        if self.results:
            self.scanpyout.write("scanpy_SCluster.h5ad",compression='gzip')
    
    # def SC_desc(self):
    #     self.desc=self.adata.copy()
    #     desc.scale(adata, zero_center=True, max_value=3)
    
    







