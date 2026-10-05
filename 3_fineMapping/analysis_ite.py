"""
gets dataframe from all information of 1MB

1.get variance for all probes
2.normalized variance
3.get related mQTLs for clumps iteration wise
"""
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

path_save="LD_fineMap/4run_regression/"

normVal=3.474051424481396
var_thr=0.05 
pip_thr=0.95 
beta_thr=0.01
ite_thr=0.05 #snp has to be selected at least 5% of the times
# %%
def identifyPatternIteration(C,M,F,I,variances):
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
            pattern="Complex" #paternal
        case [0,1]:
            pattern="Complex" #maternal
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


def pip_perCluster(ind_mqtl,rsid_mqtl,pip_thr,variances): #test ind_mqtl=410;rsid_mqtl="rs182866441";cluster=clusters[137];variances=[1,1,1,1];ite=137
    ind_mqtl=str(ind_mqtl)
    rsid_mqtl=str(rsid_mqtl)
    beta_ite_path="LD_fineMap/4run_regression/out_ite/"+ind_mqtl+"/beta_iterations.npy"
    beta_ite=np.load(beta_ite_path)

    pip_ite_path="LD_fineMap/4run_regression/out_ite/"+ind_mqtl+"/pip_iterations.npy"
    pip_ite=np.load(pip_ite_path)

    clusters_path="LD_fineMap/3.5_1Mb/clumps/"+rsid_mqtl+".txt"
    clusters=[]
    infile = open(clusters_path,'r')
    for line in infile:
        clusters.append(line.strip().split(','))
    infile.close()

    snps_path="LD_fineMap/3matrix/X_matrix/out_x/snps_"+rsid_mqtl+".txt"
    snps=np.loadtxt(snps_path,dtype=str)
    snps=list(snps[:,1])

    list_patterns_imp=["Paternal","Maternal","Complex"]
    list_patterns=["Direct","Indirect maternal","Indirect paternal"]
    dict_imp={x:[] for x in list_patterns_imp}
    dict_pat={x:[] for x in list_patterns}
    
    ind_snp=[];rsids=[];iterations=[]
    no_clusters=[];cluster_pip=[]
    beta_c=[];beta_m=[];beta_f=[];beta_i=[]
    for no_clust in range(0,len(clusters)):
        cluster=clusters[no_clust]
        ind_cluster=[snps.index(x) for x in cluster] #get indexes of snsp in cluster
        # pip per cluster
        pip_cluster=pip_ite[:,ind_cluster] #filter pip iteration for snps in cluster
        included=(np.any(pip_cluster>0,axis=1)*1).tolist() #inlusion per iteration
        calc_pip=np.sum(included)/pip_ite.shape[0]
        if calc_pip>=pip_thr: #return only clusters that passed thr
            #get pattern per iteration for each SNP
            beta_cluster=beta_ite[:,ind_cluster]
            included_snps=np.where(pip_ite[:,ind_cluster]!=0)
            no_snps=0
            #other way around, per selected snp check every iteration
            for snp in np.unique(included_snps[1]):
                entries=np.where(included_snps[1]==snp)[0].tolist()
                ite_snp=included_snps[0][entries].tolist()
                if len(ite_snp)>=int(pip_ite.shape[0]*ite_thr): #if the snp in cluster is there at least 5% of iterations of cluster
                    no_snps+=1
                    imp_pats=[];no_imp_pats=[]
                    beta_mean=beta_cluster[:,snp]
                    beta_mean[beta_mean==0]=np.nan
                    beta_mean=np.nanmean(beta_mean,axis=0)
                    for ite in ite_snp:
                        beta_snp=beta_cluster[ite,snp]
                        beta_snp=[np.sign(x) if np.abs(x)>=beta_thr else 0 for x in beta_snp] #get sign if passed beta threshold
                        imp_pat,pat=identifyPatternIteration(beta_snp[0],beta_snp[1],
                                                             beta_snp[2],beta_snp[3],variances)
                        imp_pats.append(imp_pat);no_imp_pats.append(pat)
                    # get number of each pattern 
                    no_imp_pats=[x for xs in no_imp_pats for x in xs] #flat list
                    for i in list_patterns_imp:
                        dict_imp[i].append(imp_pats.count(i))
                    for i in list_patterns:
                        dict_pat[i].append(no_imp_pats.count(i))
                    
                    #create columns for df
                    ind_snp.append(ind_cluster[snp])
                    rsids.append(cluster[snp])
                    iterations.append(len(ite_snp)) #how many times is the snp included
                    beta_c.append(beta_mean[0]);beta_m.append(beta_mean[1])
                    beta_f.append(beta_mean[2]);beta_i.append(beta_mean[3])
            #per cluster
            no_clusters+=no_snps*[no_clust] 
            cluster_pip+=no_snps*[calc_pip] #pip of cluster
             
    dict_all={}
    dict_all.update(dict_imp);dict_all.update(dict_pat)
    dict_all["Index cluster SNP"]=ind_snp;dict_all["rsID cluster SNP"]=rsids
    dict_all["Iterations"]=iterations
    dict_all["Cluster"]=no_clusters;dict_all["PIP cluster"]=cluster_pip
    dict_all["Original mQTL"]=[rsid_mqtl]*len(rsids)
    dict_all["beta C"]=beta_c;dict_all["beta M"]=beta_m;dict_all["beta F"]=beta_f
    dict_all["beta I"]=beta_i

    df=pd.DataFrame(dict_all)
    return(df)

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


def get_mQTLs(pip,beta,path_snps,var_dir):
    snps=np.genfromtxt(path_snps,dtype=object)
    ind_pip=np.where(pip>pip_thr)[0].tolist()
    filt_beta=np.array(np.where(np.abs(beta)>beta_thr))
    ind_beta=filt_beta[0,:].tolist()
    
    ind_both=list(set(ind_pip)&set(ind_beta))
    snps_name=[snps[i,1].decode() for i in ind_both]
    snps_loc=[snps[i,0].decode() for i in ind_both]
    
    #variance of coefficients
    zf = zipfile.ZipFile(var_dir)
    var_beta=np.genfromtxt(zf.open('var_beta.csv'), delimiter=',')
    var_beta=np.delete(var_beta,(0),axis=0) #intercept value in 0
    
    var_snp=[var_beta[i,:] for i in ind_both]
    var_snp=np.array(var_snp)
    #var_snp=var_beta[ind_snp,:]
    
    return(ind_both,snps_name,snps_loc,pip[ind_both].tolist(),beta[ind_both],var_snp)
    
    
path_jodie="LD_fineMap/4run_regression/out_ite/"
path_name_probes="screen_full/X_matrix/out_methylation/probes_chr_"

path_mQTLs="LD_fineMap/3.5_1Mb/mQTLs.txt"
file_mQTLs=pd.read_csv(path_mQTLs,header=None,sep="\t")
# %%
########
rows_probes=[];rows_filtered=[];noConverged_index=[]
ind_df_passed=[];w=0
cluster_df=pd.DataFrame()
#get dataframe from all probes in second screening
for i,row in file_mQTLs.iterrows():
    chro=row[1]
    if chro==23:
        print("23")
    ind=i+1
    rsid=row[0]
    ind_probe=row[3]
    folders=path_jodie+str(ind)+"/"
    
    if os.path.exists(folders+"mean_beta.csv.zip"):
        var,beta,pip,var_v=loadJODIEResults(folders)
        probe_id,loc=getLocationProbe(ind_probe,path_name_probes+str(chro)+".csv")
        #normalize variance
        norm_var=var/normVal
        
        this_probe_row=[chro,ind_probe,probe_id,loc,rsid,norm_var[0,0],norm_var[1,1],norm_var[2,2],norm_var[3,3], #normalized variance
                            var_v[0,0],var_v[1,1],var_v[2,2],var_v[3,3],norm_var[0,1],norm_var[0,2],norm_var[0,3], #variance of variances and covariances C
                           norm_var[1,2],norm_var[1,3],norm_var[2,3],var_v[0,1],var_v[0,2],var_v[0,3],var_v[1,2], #covariances and variances of covariances
                           var_v[1,3],var_v[2,3],var[0,0],var[1,1],var[2,2],var[3,3],var[0,1],var[0,2],var[0,3], #raw variances
                           var[1,2],var[1,3],var[2,3]]
        
        ###filter mQTLs
        if np.all(var_v<normVal):
            this_probe_row.append(True)
            ### get clusters and pip of clusters 
            df_mQTL=pip_perCluster(ind,rsid,pip_thr,[norm_var[0,0],norm_var[1,1],
                                             norm_var[2,2],norm_var[3,3]])
            df_mQTL["Probe ID"]=len(df_mQTL["Original mQTL"])*[probe_id]
            df_mQTL["Probe position"]=len(df_mQTL["Original mQTL"])*[loc]
            df_mQTL["Index probe"]=len(df_mQTL["Original mQTL"])*[ind_probe]
            df_mQTL["Chromosome"]=len(df_mQTL["Original mQTL"])*[chro]
            df_mQTL["VarC"]=len(df_mQTL["Original mQTL"])*[norm_var[0,0]]
            df_mQTL["VarM"]=len(df_mQTL["Original mQTL"])*[norm_var[1,1]]
            df_mQTL["VarF"]=len(df_mQTL["Original mQTL"])*[norm_var[2,2]]
            df_mQTL["VarI"]=len(df_mQTL["Original mQTL"])*[norm_var[3,3]]
            
            cluster_df=pd.concat([cluster_df,df_mQTL])#ignore_index=True
        else:
            this_probe_row.append(False)
        rows_probes.append(this_probe_row)
    else:
        noConverged_index.append(i)
        print('No converged: ',ind,chro,rsid,ind_probe)
            
noConverged=file_mQTLs[file_mQTLs.index.isin(noConverged_index)]
noConverged['newIndex']=[x+1 for x in noConverged_index]
noConverged=noConverged.set_index('newIndex')
noConverged.to_csv(path_save+'noConverged.txt',sep=' ',index=True,header=False)

#make df
col_name=["Chromosome","Index","Probe ID","Position","RSID original","varC","varM","varF","varI","var varC","var varM","var varF","var varI",
             "covCM","covCF","covCI","covMF","covMI","covFI","var covCM","var covCF","var covCI","var covMF",
             "var covMI","var covFI","raw C","raw M","raw F","raw I","raw CM",
             "raw CF","raw CI","raw MF","raw MI","raw FI","Passed variance of variances"]

probes_df=pd.DataFrame(rows_probes,columns=col_name)

##add anotations 
annot_df=pyreadr.read_r("GSM-preprocessed/EPIC_AnnotationObject_df.rds")[None]
annot_df=annot_df.rename(columns={"Name":"Probe ID"})
probes_df=pd.merge(probes_df,annot_df[["UCSC_RefGene_Name","Probe ID"]],on="Probe ID",how="inner")
probes_df=probes_df.rename(columns={'UCSC_RefGene_Name':'Gene'})
probes_df.to_csv(path_save+"allProbes_screenLD.csv",index=False)

probes_df.to_csv(path_save+"probes_clusters.csv",index=False)

## Clusters SNPs
##############
## add anotation
cluster_df=pd.merge(cluster_df,annot_df[["UCSC_RefGene_Name","Probe ID"]],on="Probe ID",how="inner")
cluster_df=cluster_df.rename(columns={'UCSC_RefGene_Name':'Gene Probe'})

##### compare with origianl mQTL
# %%
### add r^2 with original mQTL
original_mQTL=pd.read_csv("POE_mQTLs/outputs/mQTLs_patternAnnotate.csv")
### original identified patterns merge Imprinting pattern Non imprinting effects Gene Probe  Probe ID and RSID 
#"Gene Probe"
cluster_df=cluster_df.merge(original_mQTL[["Probe ID","RSID","Imprinting pattern","Non imprinting effects"]],left_on=["Probe ID","Original mQTL"],
                  right_on=["Probe ID","RSID"])

cluster_df.to_csv(path_save+"mQTL_cluster.csv",index=False)
print("Files saved")

r2_list=[]
for i,row in cluster_df.iterrows():
    original_snp=row["Original mQTL"]
    new_snp=row["rsID cluster SNP"]   
    if original_snp==new_snp:
        r2_list.append(1)
    else:
        r2_dir="LD_fineMap/3.5_1Mb/out/"+original_snp+".ld"
        r2_file=pd.read_csv(r2_dir,sep = "\s+|\t+|\s+\t+|\t+\s+")
        r2_file=r2_file[["SNP_A","SNP_B","R2"]]
        filt_r2=r2_file[((r2_file['SNP_A']==original_snp) & (r2_file['SNP_B']==new_snp)) | ((r2_file['SNP_A']==new_snp) & (r2_file['SNP_B']==original_snp))]
        if len(filt_r2)==0:
            r2_list.append(0)
        else:
            r2_list.append(filt_r2.iloc[0]["R2"])

cluster_df["r2 original SNP"]=r2_list

cluster_df.to_csv(path_save+"mQTL_clusterAnnotate.csv",index=False)


