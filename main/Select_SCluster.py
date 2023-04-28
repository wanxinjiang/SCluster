import pandas as pd 
import numpy as np
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import fowlkes_mallows_score
from SCluster import *
from BCubed import *

def Select_SCluster(adata,c_nums=None,sc3=False,soup=False,seurat=False,RaceID3=False,
                scanpy=False,cidr=False,simlr=False,sincera=False,sharp=False,results=False,alpha=0.5):
    
    results_h5ad=SClustering(adata,c_nums=None,sc3=False,soup=False,seurat=False,RaceID3=False,
                scanpy=False,cidr=False,simlr=False,sincera=False,sharp=False,results=False)
    
    methods_list=list(results_h5ad.obs.colums)

    if 'sc3' in methods_list:
        sc3=results_h5ad.obs['sc3']
    
    if 'soup' in methods_list:
        soup=results_h5ad.obs['soup']
    
    if 'seurat' in methods_list:
        seurat=results_h5ad.obs['seurat']
    
    if 'RaceID3' in methods_list:
        RaceID3=results_h5ad.obs['RaceID3']
    
    if 'scanpy' in methods_list:
        scanpy=results_h5ad.obs['scanpy']
    
    if 'cidr' in methods_list:
        cidr=results_h5ad.obs['cidr']
    
    if 'simlr' in methods_list:
        simlr=results_h5ad.obs['simlr']
    
    if 'sincera' in methods_list:
        sincera=results_h5ad.obs['sincera']
    
    if 'sharp' in methods_list:
        sharp=results_h5ad.obs['sharp']
    
    SC_cluster=[sc3,soup,seurat,RaceID3,scanpy,cidr,simlr,sincera,sharp]

    SC_eva={}
    for a,sc_1 in enumerate(SC_cluster): 
        for b,sc_2 in enumerate(SC_cluster[a+1:]):
            ARI=adjusted_rand_score(sc_1,sc_2)
            NMI=normalized_mutual_info_score(sc_1,sc_2)
            FM=fowlkes_mallows_score(sc_1,sc_2)
            BCuded_f1=bcubed(sc_1,sc_2)
            SC_score=alpha*(ARI+NMI+FM)+(1-alpha)*BCuded_f1
            SC_cluster[str(a)+'_'+str(b)]=SC_score
        if a==len(SC_cluster)-1:
            break
        


