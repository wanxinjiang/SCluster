import pandas as pd 
import numpy as np
import anndata
import re
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import fowlkes_mallows_score
from sklearn.metrics.cluster import homogeneity_completeness_v_measure
from sklearn.metrics.cluster import davies_bouldin_score,calinski_harabasz_score
from SClusterclass import *
from BCubed import *
from Jaccard import *
import anndata2ri
import time

def get_cluster(adata,c_nums=None,sc3=False,soup=False,seurat=False,RaceID3=False,
                scanpy=False,cidr=False,simlr=False,sincera=False,sharp=False,results=False,alpha=0.5):
    
    results_h5ad=SClustering(adata=adata,c_nums=c_nums,sc3=sc3,soup=soup,
                              seurat=seurat,RaceID3=RaceID3,scanpy=scanpy,
                              cidr=cidr,simlr=simlr,sincera=sincera,
                            sharp=sharp,results=results)
    
    methods_list=list(results_h5ad.obs.columns)
    
    _SC_cluster=[]
    if 'sc3' in methods_list:
        sc3=results_h5ad.obs['sc3']
        _SC_cluster.append(['sc3',sc3])

    if 'soup' in methods_list:
        soup=results_h5ad.obs['soup']
        _SC_cluster.append(['soup',soup])
    
    if 'seurat' in methods_list:
        seurat=results_h5ad.obs['seurat']
        _SC_cluster.append(['seurat',seurat])
    
    if 'RaceID3' in methods_list:
        RaceID3=results_h5ad.obs['RaceID3']
        _SC_cluster.append(['RaceID3',RaceID3])
    
    if 'scanpy' in methods_list:
        scanpy=results_h5ad.obs['scanpy']
        _SC_cluster.append(['scanpy',scanpy])
    
    if 'cidr' in methods_list:
        cidr=results_h5ad.obs['cidr']
        _SC_cluster.append(['cidr',cidr])
    
    if 'simlr' in methods_list:
        simlr=results_h5ad.obs['simlr']
        _SC_cluster.append(['simlr',simlr])
    
    if 'sincera' in methods_list:
        sincera=results_h5ad.obs['sincera']
        _SC_cluster.append(['sincera',sincera])
    
    if 'sharp' in methods_list:
        sharp=results_h5ad.obs['sharp']
        _SC_cluster.append(['sharp',sharp])
    
    return results_h5ad,_SC_cluster


 
def SCluster(adata,c_nums=None,\
            sc3=False,soup=False,seurat=False,RaceID3=False,\
            scanpy=False,cidr=False,simlr=False,sincera=False,sharp=False,\
            results=False,alpha=0.5,first_result=False):
    print('Performing SCluster clustering... ')
    start_time=time.time()
    
    adada_res,_SC_cluster=get_cluster(adata=adata,c_nums=c_nums,sc3=sc3,soup=soup,
                                                           seurat=seurat,RaceID3=RaceID3,scanpy=scanpy,
                                                           cidr=cidr,simlr=simlr,sincera=sincera,
                                                           sharp=sharp,results=results,alpha=alpha) 
    
    
    
    
    SCluster_results=sc_results.obs[similarity_w]
    
    if sc3:
        del sc_results.obs['sc3']
    if soup:
       del sc_results.obs['soup']
    if cidr:
        del sc_results.obs['cidr']
    if simlr:
        del sc_results.obs['simlr']
    if sincera:
        del sc_results.obs['sincera']
    if sharp:
        del sc_results.obs['sharp']
    if seurat:
        del sc_results.obs['seurat']
    if RaceID3:
        del sc_results.obs['RaceID3'] 
    if scanpy:
        del sc_results.obs['scanpy']
    
    sc_results.obs['SCluster']=SCluster_results
    
    # sc_results.obs.index=pd.DataFrame(sc_results.obs.index)[0].apply(lambda x : re.findall(r'(.*?)-.*',x)[0])
    # sc_results.obs.index.name='' 
    end_time=time.time()
    print('SCluster done...')
    print('SCluster takes {} secs'.format(end_time-start_time))
    if first_result:   
        return sc_results,results_h5ad_first
    else:
        return sc_results

    
    
    




    




    
        


