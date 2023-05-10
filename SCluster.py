import anndata
import re
import time
import pandas as pd 
import numpy as np
import scanpy as sc
import anndata2ri
import itertools
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import fowlkes_mallows_score
from sklearn.metrics.cluster import homogeneity_completeness_v_measure
from sklearn.metrics.cluster import davies_bouldin_score,calinski_harabasz_score
from .SClusterclass import *
from .BCubed import *
from .Jaccard import *
from .ClusterEnsembles import cluster_ensembles


def get_cluster(adata,c_nums=None,sc3=False,soup=False,seurat=False,RaceID3=False,
                scanpy=False,cidr=False,simlr=False,sincera=False,sharp=False,results=False,alpha=0.82):
    
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
    
    similarity_w,_SC_cluster_dict,Small_similarity=similarity_method(_SC_cluster)

    cat_adata,second_data,c_nums=adata_cat(similarity_w,_SC_cluster_dict,results_h5ad,alpha)
    
    print(similarity_w,c_nums,Small_similarity)
    
    return cat_adata,second_data,results_h5ad,similarity_w,c_nums,_SC_cluster,Small_similarity


def similarity_method(_SC_cluster):
    
    SC_eva=np.zeros((len(_SC_cluster),len(_SC_cluster)))
    _SC_flag=[i[0] for i in _SC_cluster]
    for a,sc_1 in enumerate(_SC_cluster): 
        for b,sc_2 in enumerate(_SC_cluster):
            ARI=adjusted_rand_score(sc_1[1],sc_2[1])
            NMI=normalized_mutual_info_score(sc_1[1],sc_2[1])
            FM=fowlkes_mallows_score(sc_1[1],sc_2[1])
            f1=bcubed(sc_1[1],sc_2[1])
            hcv=np.mean(homogeneity_completeness_v_measure(sc_1[1],sc_2[1]))
            SC_score=0.5*(ARI+NMI+FM)+0.5*(hcv+f1)
            SC_eva[a][b]=SC_score
    SC_eva=pd.DataFrame(SC_eva)
    SC_eva.index=_SC_flag
    SC_eva.columns=_SC_flag
    eva_mean=dict(SC_eva.mean())
    similarity_w=list(eva_mean.keys())[list(eva_mean.values()).index(max(eva_mean.values()))]
    Small_similarity=''
    if len(_SC_cluster)>=2:
        Small_similarity=list(eva_mean.keys())[list(eva_mean.values()).index(min(eva_mean.values()))]
    _SC_cluster_dict=dict(_SC_cluster)

    return similarity_w,_SC_cluster_dict,Small_similarity

def adata_cat(similarity_w,_SC_cluster_dict,adata,alpha):
    ref_way=_SC_cluster_dict[similarity_w]
    ref_cluster=list(np.unique(ref_way))
    ref_cluster_set=ref_cluster.copy()
    similarity_cluster={}
    for key,values in _SC_cluster_dict.items():
        _s=pd.DataFrame(jac(values,ref_way))
        _s=_s[np.max(_s,axis=1)>=alpha]
        _s=list(_s.index)
        similarity_cluster[key]=list(set(ref_cluster_set)&set(_s))

    cluster_sim=[i for i in similarity_cluster.values()]
    cluster_sim=list(itertools.chain.from_iterable(cluster_sim))  
    cluster_sim=pd.value_counts(cluster_sim)
    if len(_SC_cluster_dict)<=3:
        cluster_sim=cluster_sim[cluster_sim>=2].index.tolist()
    elif len(_SC_cluster_dict)<=6:
        cluster_sim=cluster_sim[cluster_sim>=abs(len(_SC_cluster_dict)-2)].index.tolist()
    else:
        cluster_sim=cluster_sim[cluster_sim>=abs(len(_SC_cluster_dict)-4)].index.tolist()
        
    cat_adata=[adata[adata.obs[similarity_w]==i] for i in cluster_sim]
    if len(cat_adata)>1:
        cat_adata=cat_adata[0].concatenate(cat_adata[1:])
        cat_adata.obs[similarity_w]=cat_adata.obs['batch'].astype('int')
    else:
        try:
            cat_adata=cat_adata[0]
            if 'batch' in cat_adata.obs.columns:
                cat_adata.obs[similarity_w]=cat_adata.obs['batch'].astype('int')
            else:
                cat_adata.obs[similarity_w]=0
        except:
            cat_adata='' 
            
    second_data=[adata[adata.obs[similarity_w]==i] for i in list(set(ref_cluster)-set(cluster_sim))]
    if len(second_data)>1:
        second_data=second_data[0].concatenate(second_data[1:])
    else:
        try:
            second_data=second_data[0]
        except:
            second_data=''
            
    c_nums=len(list(set(ref_cluster)-set(cluster_sim)))

    return cat_adata,second_data,c_nums


def Ensemble(_SC_cluster,Small_similarity):
    
    Ensemble_cluster=[np.array(i[1].astype('int') ) for i in _SC_cluster if i[0]!=Small_similarity ]
    Ensemble_cluster=cluster_ensembles(np.array(Ensemble_cluster),solver='all')
    
    return Ensemble_cluster
    

def SCluster(adata,c_nums=None,\
            sc3=False,soup=False,seurat=False,RaceID3=False,\
            scanpy=False,cidr=False,simlr=False,sincera=False,sharp=False,\
            results=False,alpha=0.82):
    
    print('Performing SCluster clustering... ')
    input_c=c_nums
    adata.obs['sort']=np.array([i for i in range(adata.n_obs)])
    start_time=time.time()
    cat_resluts,second_data,results_h5ad_first,\
    similarity_w,c_nums,_SC_cluster,Small_similarity=get_cluster(adata=adata,c_nums=c_nums,\
                                                            sc3=sc3,soup=soup,\
                                                           seurat=seurat,RaceID3=RaceID3,scanpy=scanpy,\
                                                           cidr=cidr,simlr=simlr,sincera=sincera,\
                                                           sharp=sharp,results=results,alpha=alpha)
    
    print(second_data)
    print('*'*50)
    print(cat_resluts) 
    print(Small_similarity)
    if Small_similarity=='sc3':
        sc3=False

    if Small_similarity=='seurat':
        seurat=False

    if Small_similarity=='RaceID3' or 'RaceID3' not in results_h5ad_first.obs.columns:
        RaceID3=False

    if Small_similarity=='soup':
        soup=False

    if Small_similarity=='scanpy':
        scanpy=False

    if Small_similarity=='cidr':
        cidr=False

    if Small_similarity=='sincera':
        sincera=False

    if Small_similarity=='simlr' or 'simlr' not in results_h5ad_first.obs.columns:
        simlr=False

    if Small_similarity=='sharp':
        sharp=False
      
    sc_results=cat_resluts
     
    try: 
        cat_nums=len(np.unique(sc_results.obs[similarity_w]))
        if input_c==c_nums:
            Ensemble_cluster=Ensemble(_SC_cluster,Small_similarity)
            sc_results=results_h5ad_first.copy()
            sc_results.obs[similarity_w]=Ensemble_cluster.values
        elif second_data.n_obs>50:
            cat_resluts,second_data,results_h5ad,similarity_s,\
            c_nums,_SC_cluster,Small_similarity=get_cluster(adata=second_data,\
                                                             c_nums=c_nums,sc3=sc3,soup=soup,\
                                                            seurat=seurat,RaceID3=RaceID3,scanpy=scanpy,\
                                                            cidr=cidr,simlr=simlr,sincera=sincera,\
                                                            sharp=sharp,results=results,alpha=alpha)
            
            
            Ensemble_cluster=Ensemble(_SC_cluster,Small_similarity)
            # Ensemble_cluster=results_h5ad.obs[similarity_s]

            results_h5ad.obs[similarity_w]=Ensemble_cluster.astype('int')+int(cat_nums)
            sc_results=sc_results.concatenate(results_h5ad)
        else:
            Ensemble_cluster=Ensemble(_SC_cluster,Small_similarity)
            # Ensemble_cluster=results_h5ad.obs[similarity_w]
            sc_results=results_h5ad_first.copy()
            sc_results.obs[similarity_w]=Ensemble_cluster.values
    except:
        Ensemble_cluster=results_h5ad_first.obs[similarity_w]
        sc_results=results_h5ad_first.copy()
        sc_results.obs[similarity_w]=Ensemble_cluster.values
     
    SCluster_results=sc_results.obs[similarity_w]

    sc_results.obs['SCluster']=SCluster_results
    sc_results.obs['SCluster']=sc_results.obs['SCluster'].astype('category')
    
    sc_results.obs = sc_results.obs.sort_values('sort', ascending=True)
    
    results_h5ad_first.obs['SCluster']=sc_results.obs['SCluster'].values
        
    if 'sort' in sc_results.obs.columns:
        del results_h5ad_first.obs['sort']
        
    sc.pp.neighbors(results_h5ad_first, n_neighbors=10, n_pcs=40)
    sc.tl.umap(results_h5ad_first)
    
    end_time=time.time()
    print('SCluster done...')
    print('SCluster takes {} secs'.format(end_time-start_time))
    
    return results_h5ad_first

    
    
    




    




    
        


