"""
gets dataframe from all information of second screening CMFI model

1.get variance for all probes
2.normalized variance
3.get related mQTLs
3.save df
4.save information probes that weren't calculated
"""

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

path_save="POE_mQTLs/outputs/"

normVal=3.474051424481396
var_thr=0.05 #varI
pip_thr=0.95 #pip threshold
beta_thr=0.01

def loadJODIEResults(JODIE_path):
    var=np.loadtxt(JODIE_path+"mean_V.txt")
    try:
      zf = zipfile.ZipFile(JODIE_path+"mean_beta.csv.zip")
    except:
      print(JODIE_path)
    beta_j=np.genfromtxt(zf.open('mean_beta.csv'), delimiter=',')
    beta_j=np.delete(beta_j,(0),axis=0) #intercept value in 0
    var_v=np.loadtxt(JODIE_path+"var_V.txt")
    
    pip_j=np.loadtxt(JODIE_path+"mean_prob.txt") #PER SNP NOT PER VARIABLE 
    return(var,beta_j,pip_j,var_v)

def getLocationProbe(probe,path):
    probes=np.genfromtxt(path,delimiter=",",dtype=object)
    
    name=probes[probe,0].decode()
    loc=probes[probe,1].decode()
    return(name,loc)

def get_mQTLs(pip,beta,path_snps):
    snps=np.genfromtxt(path_snps,dtype=object)
    ind_pip=np.where(pip>pip_thr)[0].tolist()
    filt_beta=np.array(np.where(np.abs(beta)>beta_thr))
    ind_beta=filt_beta[0,:].tolist()
    
    ind_both=list(set(ind_pip)&set(ind_beta))
    snps_name=[snps[i,1].decode() for i in ind_both]
    snps_loc=[snps[i,0].decode() for i in ind_both]
    return(ind_both,snps_name,snps_loc,pip[ind_both].tolist(),beta[ind_both])
    
    
path_jodie="screen_full/"
path_name_probes="screen_full/X_matrix/out_methylation/probes_chr_"
path_snps="screen_full/X_matrix/out_x/snps_"

############################################
rows_probes=[];rows_mQTLs=[];noConverged=[]
ind_df_passed=[];w=0
#get dataframe from all probes in second screening
for chro_int in range(0,24):
    chro=str(chro_int)
    folders=os.walk(path_jodie+chro+"/jodie/out/")
    for idx,(dirs,null,files) in enumerate(folders):
        if idx==0: continue #first iteration gives main path
        ind_probe=int(dirs.split("/")[len(dirs.split("/"))-1])
        #make sure probe converged
        if "mean_beta.csv.zip" in files:
            var,beta,pip,var_v=loadJODIEResults(dirs+"/")
            probe_id,loc=getLocationProbe(ind_probe,path_name_probes+chro+".csv")
            #normalize variance
            norm_var=var/normVal
            
            rows_probes.append([chro,ind_probe,probe_id,loc,norm_var[0,0],norm_var[1,1],norm_var[2,2],norm_var[3,3], #normalized variance
                                var_v[0,0],var_v[1,1],var_v[2,2],var_v[3,3],norm_var[0,1],norm_var[0,2],norm_var[0,3], #variance of variances and covariances C
                               norm_var[1,2],norm_var[1,3],norm_var[2,3],var_v[0,1],var_v[0,2],var_v[0,3],var_v[1,2], #covariances and variances of covariances
                               var_v[1,3],var_v[2,3],var[0,0],var[1,1],var[2,2],var[3,3],var[0,1],var[0,2],var[0,3], #raw variances
                               var[1,2],var[1,3],var[2,3]])
            
            ###filter mQTLs
            if norm_var[3,3]>var_thr and np.all(var_v<normVal): #varI>0.05 and variances of variances < norm value
                ind_df_passed.append(w)#index in df for probes that passed filters
                snps_ind,snps_name,snps_loc,snps_pip,snps_beta=get_mQTLs(pip, beta, path_snps+chro+".txt")
                
                if len(snps_ind)>0: #there are mQTLs in this probe
                    rows_thisProbe=[[chro]*len(snps_ind),[ind_probe]*len(snps_ind),[probe_id]*len(snps_ind),[loc]*len(snps_ind),
                                    [norm_var[0,0]]*len(snps_ind),[norm_var[3,3]]*len(snps_ind),[norm_var[0,3]]*len(snps_ind),
                                    snps_ind,snps_name,snps_loc,snps_pip,snps_beta[:,0].tolist(),snps_beta[:,1].tolist(),
                                    snps_beta[:,2].tolist(),snps_beta[:,3].tolist()]
                    if len(rows_mQTLs)==0:#first entry
                        rows_mQTLs=rows_thisProbe
                    else:
                        rows_mQTLs=[a+b for a,b in zip(rows_mQTLs,rows_thisProbe)]#append lists
            w+=1
        else:
            noConverged.append([chro,ind_probe])
#make df
col_name=["Chromosome","Index","Probe ID","Position","varC","varM","varF","varI","var varC","var varM","var varF","var varI",
             "covCM","covCF","covCI","covMF","covMI","covFI","var covCM","var covCF","var covCI","var covMF",
             "var covMI","var covFI","raw C","raw M","raw F","raw I","raw CM",
             "raw CF","raw CI","raw MF","raw MI","raw FI"]

allProbesSecond=pd.DataFrame(rows_probes,columns=col_name)
##add anotations 
annot_df=pyreadr.read_r("EPIC_AnnotationObject_df.rds")[None]
annot_df=annot_df.rename(columns={"Name":"Probe ID"})
allProbesSecond=pd.merge(allProbesSecond,annot_df[["UCSC_RefGene_Name","Probe ID"]],on="Probe ID",how="inner")
allProbesSecond=allProbesSecond.rename(columns={'UCSC_RefGene_Name':'Gene'})
allProbesSecond.to_csv(path_save+"screen2CMFI.csv",index=False)


passedFilters=allProbesSecond[allProbesSecond.index.isin(ind_df_passed)]
passedFilters.to_csv(path_save+"poe_probes.csv",index=False)

#SNPs
col_mQTLs=["Chromosome","Index probe","Probe ID","Position Probe","varC","varI","covCI","Index SNP","RSID","Position SNP","PIP","beta C",
           "beta M","beta F","beta I"]
mQTLs=pd.DataFrame(rows_mQTLs)
mQTLs=mQTLs.transpose()
mQTLs.columns=col_mQTLs
mQTLs.to_csv(path_save+"mQTLs.csv",index=False)
print("Files saved")

