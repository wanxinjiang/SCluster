import anndata
import re
import time
import pandas as pd 
import numpy as np
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import fowlkes_mallows_score
from sklearn.metrics.cluster import homogeneity_completeness_v_measure
from sklearn.metrics.cluster import davies_bouldin_score,calinski_harabasz_score
from SClusterclass import *
from BCubed import *
from Jaccard import *
import anndata2ri
import itertools
import ClusterEnsembles as CE

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
    
    similarity_w,_SC_cluster_dict,Small_similarity=similarity_method(_SC_cluster,alpha)

    cat_adata,second_data,c_nums=adata_cat(similarity_w,_SC_cluster_dict,results_h5ad)
    
    print(similarity_w,c_nums,Small_similarity)
    
    return cat_adata,second_data,results_h5ad,similarity_w,c_nums,_SC_cluster,Small_similarity


def similarity_method(_SC_cluster,alpha):
    
    SC_eva=np.zeros((len(_SC_cluster),len(_SC_cluster)))
    _SC_flag=[i[0] for i in _SC_cluster]
    for a,sc_1 in enumerate(_SC_cluster): 
        for b,sc_2 in enumerate(_SC_cluster):
            ARI=adjusted_rand_score(sc_1[1],sc_2[1])
            NMI=normalized_mutual_info_score(sc_1[1],sc_2[1])
            FM=fowlkes_mallows_score(sc_1[1],sc_2[1])
            f1=bcubed(sc_1[1],sc_2[1])
            hcv=np.mean(homogeneity_completeness_v_measure(sc_1[1],sc_2[1]))
            SC_score=alpha*(ARI+NMI+FM)+(1-alpha)*(hcv+f1)
            SC_eva[a][b]=SC_score
    SC_eva=pd.DataFrame(SC_eva)
    SC_eva.index=_SC_flag
    SC_eva.columns=_SC_flag
    eva_mean=dict(SC_eva.mean())
    similarity_w=list(eva_mean.keys())[list(eva_mean.values()).index(max(eva_mean.values()))]
    if len(_SC_cluster)>=2:
        Small_similarity=list(eva_mean.keys())[list(eva_mean.values()).index(min(eva_mean.values()))]
    _SC_cluster_dict=dict(_SC_cluster)

    return similarity_w,_SC_cluster_dict,Small_similarity

def adata_cat(similarity_w,_SC_cluster_dict,adata):
    ref_way=_SC_cluster_dict[similarity_w]
    ref_cluster=list(np.unique(ref_way))
    ref_cluster_set=ref_cluster.copy()
    similarity_cluster={}
    for key,values in _SC_cluster_dict.items():
        _s=pd.DataFrame(jac(values,ref_way))
        _s=_s[np.max(_s,axis=1)>=0.75]
        _s=list(_s.index)
        similarity_cluster[key]=list(set(ref_cluster_set)&set(_s))

    cluster_sim=[i for i in similarity_cluster.values()]
    cluster_sim=list(itertools.chain.from_iterable(cluster_sim))  
    cluster_sim=pd.value_counts(cluster_sim)
    cluster_sim=cluster_sim[cluster_sim>=len(_SC_cluster_dict)-2].index.tolist()
        
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


def get_subcluster(sub_results,_SC_cluster):
    try:
        X=sub_results.X.toarray()
    except:
        X=sub_results.X
    dc_values={}
    for i in range(len(_SC_cluster)):
        db_score=davies_bouldin_score(X,_SC_cluster[i][1])
        ch_score=calinski_harabasz_score(X,_SC_cluster[i][1])
        scores=np.log10(db_score)+np.log10(ch_score)
        dc_values[_SC_cluster[i][0]]=scores
    sub_w=list(dc_values.keys())[list(dc_values.values()).index(max(dc_values.values()))]
    return sub_results.obs[sub_w]

def Ensemble(_SC_cluster,Small_similarity):
    
    Ensemble_cluster=[np.array(i[1].astype('int') ) for i in _SC_cluster if i[0]!=Small_similarity ]
    Ensemble_cluster=CE.cluster_ensembles(np.array(Ensemble_cluster),solver='mcla')
    
    return Ensemble_cluster
    

def SCluster(adata,c_nums=None,\
            sc3=False,soup=False,seurat=False,RaceID3=False,\
            scanpy=False,cidr=False,simlr=False,sincera=False,sharp=False,\
            results=False,alpha=0.5,first_result=False):
    
    print('Performing SCluster clustering... ')
    input_c=c_nums
    start_time=time.time()
    cat_resluts,second_data,results_h5ad_first,\
    similarity_w,c_nums,_SC_cluster,Small_similarity=get_cluster(adata=adata,c_nums=c_nums,\
                                                            sc3=sc3,soup=soup,\
                                                           seurat=seurat,RaceID3=RaceID3,scanpy=scanpy,\
                                                           cidr=cidr,simlr=simlr,sincera=sincera,\
                                                           sharp=sharp,results=results,alpha=alpha)
    
    print(second_data)
    print('bug0')
    print('*'*50)
    
    sc3_flag=True
    seurat_flag=True
    RaceID3_flag=True
    soup_flag=True
    scanpy_flag=True
    cidr_flag=True
    sincera_flag=True
    simlr_flag=True
    sharp_flag=True
    
    if Small_similarity=='sc3':
        sc3=False
        sc3_flag=False
    
    if Small_similarity=='seurat':
        seurat=False
        seurat_flag=False
    
    if Small_similarity=='RaceID3':
        RaceID3=False
        RaceID3_flag=False
    
    if Small_similarity=='soup':
        soup=False
        soup_flag=False
    
    if Small_similarity=='scanpy':
        scanpy=False
        scanpy_flag=False
    if Small_similarity=='cidr':
        cidr=False
        cidr_flag=False
    
    if Small_similarity=='sincera':
        sincera=False
        sincera_flag=False
    
    if Small_similarity=='simlr':
        simlr=False
        simlr_flag=False
    
    if Small_similarity=='sharp':
        sharp=False
        sharp_flag=False
        
    sc_results=cat_resluts
    print(sc_results)  
    try: 
        cat_nums=len(np.unique(sc_results.obs[similarity_w]))
        if input_c==c_nums:
            Ensemble_cluster=Ensemble(_SC_cluster,Small_similarity)
            sc_results=results_h5ad_first.copy()
            sc_results.obs[similarity_w]=Ensemble_cluster.values
        if second_data.n_obs>50:
            cat_resluts,second_data,results_h5ad,similarity_s,\
            c_nums,_SC_cluster,Small_similarity=get_cluster(adata=second_data,\
                                                             c_nums=c_nums,sc3=sc3,soup=soup,\
                                                            seurat=seurat,RaceID3=RaceID3,scanpy=scanpy,\
                                                            cidr=cidr,simlr=simlr,sincera=sincera,\
                                                            sharp=sharp,results=results,alpha=alpha)
            
            print('bug1')
            Ensemble_cluster=Ensemble(_SC_cluster,Small_similarity)
            # Ensemble_cluster=results_h5ad.obs[similarity_s]
            print('bug2')
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
    
    test=''
    try:
        test=results_h5ad
    except:
        pass
    
    
    if sc3 and sc3_flag:
        del sc_results.obs['sc3']
    if soup and soup_flag:
       del sc_results.obs['soup']
    if cidr and cidr_flag:
        del sc_results.obs['cidr']
    if simlr and simlr_flag:
        del sc_results.obs['simlr']
    if sincera and sincera_flag:
        del sc_results.obs['sincera']
    if sharp and sharp_flag:
        del sc_results.obs['sharp']
    if seurat and seurat_flag:
        del sc_results.obs['seurat']
    if RaceID3 and RaceID3_flag:
        del sc_results.obs['RaceID3'] 
    if scanpy and scanpy_flag:
        del sc_results.obs['scanpy']
    
    sc_results.obs['SCluster']=SCluster_results
    
    # sc_results.obs.index=pd.DataFrame(sc_results.obs.index)[0].apply(lambda x : re.findall(r'(.*?)-.*',x)[0])
    # sc_results.obs.index.name='' 
    end_time=time.time()
    print('SCluster done...')
    print('SCluster takes {} secs'.format(end_time-start_time))
    if first_result:   
        return sc_results,results_h5ad_first,test
    else:
        return sc_results

    
    
    




    




    
        


