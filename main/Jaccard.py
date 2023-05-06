import json
import sys
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial import distance
import pandas as pd
import numpy as np
def jac(method,df):
    df=pd.DataFrame(df)
    df.columns=['type']
    method=pd.DataFrame(method)
    method.columns=['cell']
    method.index=df.index
    dict1={}
    for i in np.unique(method['cell']):
        dict1[i]=method[method["cell"]==i].index.tolist()
    dict2={}
    for i in np.unique(df['type']):
        dict2[i]=df[df["type"]==i].index.tolist()
        
    row = list(dict1.keys())
    col = list(dict2.keys())
    row.sort()
    col.sort()
    sims = np.zeros((len(row),len(col)))
    def intersection_list(list1, list2):  
        list3 = [value for value in list1 if value in list2]  
        return len(list3)
    def union_list(list1, list2):
        mm ={}
        for x in list1:
            mm[x]=1
        for y in list2:
            mm[y]=1
        return len(mm)

    for ri,rname in enumerate(row):
        for ci,cname in enumerate(col):
            genes1 = dict1[rname]
            genes2 = dict2[cname]
            if len(genes1)>0 and len(genes2)>0:
                sims[ri,ci] = intersection_list(genes1,genes2)/union_list(genes1,genes2)
    ret = pd.DataFrame(data=sims,columns=col,index=row)
    ret['max']=ret.max(axis=1)
    ret=ret.sort_values(by=['max'],ascending=False)
    del ret['max']
    ret=ret.T
    ret['max']=ret.max(axis=1).values
    ret=ret.sort_values(by=['max'],ascending=False)
    del ret['max']
    return ret