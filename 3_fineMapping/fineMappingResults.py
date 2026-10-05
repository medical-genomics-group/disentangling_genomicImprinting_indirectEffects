# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import zarr
import os
import sys
import glob
import zipfile
from itertools import compress
import csv
import json
import pyreadr
import pickle

var_thr=0.05
beta_thr=0.01

def identifyPattern(C,M,F,I,variances):
    if variances[0]<var_thr:
        C=0
    if variances[1]<var_thr:
        M=0
    if variances[2]<var_thr:
        F=0
    if variances[3]<var_thr:
        I=0
        
    signs=[C,I]
    match signs:
        #imprinting
        case [-1,-1]:
            pattern="Paternal" 
        case [1,-1,]:
            pattern="Maternal" 
        case [-1,1]:
            pattern="Maternal" 
        case [1,1]:
            pattern="Paternal" 
        case [0,-1]:
            pattern="Complex" 
        case [0,1]:
            pattern="Complex"
        case _:
            pattern=np.nan
   
    other=[] 
    if np.abs(C)!=0:
        other.append("Direct")
    if np.abs(M)!=0:
        other.append("Indirect maternal")
    if np.abs(F)!=0:
        other.append("Indirect paternal")
    
    return(pattern,other)

def signCoeff(coeff):
    if np.abs(coeff)>=beta_thr:
        return(np.sign(coeff))
    else:
        return(0)

def signEffectsSNP(df,cols):
    for col in cols:
        df["sign_"+col]=df.apply(lambda row: signCoeff(["beta "+col]),axis=1)

clusters_df_path="LD_fineMap/4run_regression/mQTL_clusterAnnotate.csv"
clusters_df=pd.read_csv(clusters_df_path)
print("How many originally, wihtou varI filter: ",len(clusters_df))

clusters_df['r2 original SNP']=clusters_df['r2 original SNP'].replace(0, np.nan)
clusters_df=clusters_df.drop(columns=['RSID'])

clusters_df=clusters_df[clusters_df["VarI"]>=0.05]
print("How many have varI>=0.05: ", len(clusters_df))

#%%
### statistics
freq_clusters=clusters_df.groupby(["Original mQTL","Index cluster SNP"]).size()
print("Number of SNPs slected by clusters: ",len(clusters_df))
print("Number of unique clusters",len(freq_clusters))
print("Number of unique original mQTL with at least one SNP selected in clusters: ",len(np.unique(clusters_df["Original mQTL"])))
print("Number of unique probes with at least one SNP selected in clusters: ",len(np.unique(clusters_df["Probe ID"])))
print("Number of original snps in clusters: ",len(clusters_df[clusters_df["Same SNP"]==True]))

### original mQTLs
mQTLs=pd.read_csv("POE_mQTLs/outputs/mQTLs_patternAnnotate.csv")

# number of things that did not run
missing=[]
for i,row in mQTLs.iterrows():
    rsid=row["RSID"]
    probe=row["Probe ID"]
    
    in_clusters=clusters_df[(clusters_df["Original mQTL"]==rsid) & (clusters_df["Probe ID"]==probe)]
    if len(in_clusters)==0:
        missing.append([rsid,probe])
        
print("Out of the 412 mQTLs how many do we get: ",412-len(missing))

###########################
#imprinting patterns
pats=["Maternal","Paternal","Complex"]
selected_pat=[]
for i,row in clusters_df.iterrows():
    iter_row=row["Iterations"]
    patterns=[row[pats[0]],row[pats[1]],row[pats[2]]]
    sign_pat=[pats[i] for i in range(0,3) if patterns[i]>=(iter_row*0.75)]
    if len(sign_pat)>0:
        selected_pat.append(sign_pat[0])
    else:
        selected_pat.append("Undefined")
        
clusters_df["Selected pattern"]=selected_pat
clusters_df["Match patterns"]=clusters_df["Selected pattern"]==clusters_df["Imprinting pattern"]

path_save="POE_mQTLs/outputs/"

#### same snp-probe association from different windows collapsed
mQTLs=clusters_df
once=mQTLs[~mQTLs.duplicated(subset=["rsID cluster SNP","Probe ID"],keep=False)]
dup=mQTLs[mQTLs.duplicated(subset=["rsID cluster SNP","Probe ID"],keep=False)]
names_dup=dup[["rsID cluster SNP","Probe ID"]].drop_duplicates()

pats=["Maternal","Paternal","Complex"]

for i,row in names_dup.iterrows():
    snp=row["rsID cluster SNP"]
    probe=row["Probe ID"]
    
    same=mQTLs[(mQTLs["rsID cluster SNP"]==snp)&(mQTLs["Probe ID"]==probe)]
    sum_same=same.sum(numeric_only=True)
    add=sum_same.copy()
    add=add.to_frame().T
    
    sign=[pats[i] for i in range(0,3) if sum_same[pats[i]]>=(sum_same["Iterations"]*0.75)]
    if len(sign)>0:
        add["Selected pattern"]=(sign[0])
    else:
        add["Selected pattern"]="Undefined"
    
    add["rsID cluster SNP"]=snp; add["Probe ID"]=probe; 
    add["Original mQTL"]=";".join(same["Original mQTL"].tolist())
    add["Gene Probe"]=same.iloc[0]["Gene Probe"]
    add["Probe position"]=same.iloc[0]["Probe position"]
    add["Index probe"]=same.iloc[0]["Index probe"]
    add["Chromosome"]=same.iloc[0]["Chromosome"]
    add["Index cluster SNP"]=";".join([str(int(x)) for x in same["Index cluster SNP"]])
    add["Cluster"]=";".join([str(int(x)) for x in same["Cluster"]])
    add["r2 original SNP"]=";".join([str(x) for x in same["r2 original SNP"]])
    add["Imprinting pattern"]=same.iloc[0]["Imprinting pattern"]
    add["PIP cluster"]=";".join([str(x) for x in same["PIP cluster"]])
    
    #weighted mean
    add["beta C"]=np.sum(same["beta C"]*same["Iterations"])/add["Iterations"]
    add["beta M"]=np.sum(same["beta M"]*same["Iterations"])/add["Iterations"]
    add["beta F"]=np.sum(same["beta F"]*same["Iterations"])/add["Iterations"]
    add["beta I"]=np.sum(same["beta I"]*same["Iterations"])/add["Iterations"]
    
    add["VarC"]=np.sum(same["VarC"]*same["Iterations"])/add["Iterations"]
    add["VarM"]=np.sum(same["VarM"]*same["Iterations"])/add["Iterations"]
    add["VarF"]=np.sum(same["VarF"]*same["Iterations"])/add["Iterations"]
    add["VarI"]=np.sum(same["VarI"]*same["Iterations"])/add["Iterations"]
    
    once=pd.concat([once,add])
    

once.to_csv(path_save+"cluster_SNPs.csv",index=False)
clusters_df.to_csv(path_save+"cluster_df.csv",index=False)

