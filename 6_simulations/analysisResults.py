"""
load results from simulations,

Get ROC values for JODIE, SuSIE and LASSO

create summary dataframe and methods in pickle format

"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import zarr
import os
import sys
import glob
import zipfile
from itertools import groupby
import csv
import pickle
import seaborn as sns

noConverge=[]

def loadJODIEResults(JODIE_path,beta_thr=False,pip_thr=False,var_thr=False,var_comp=True):
    try:
        covars=np.loadtxt(JODIE_path+"mean_V.txt")
        zf = zipfile.ZipFile(JODIE_path+"mean_beta.csv.zip")
        beta_j=np.genfromtxt(zf.open('mean_beta.csv'), delimiter=',')
        beta_j=np.delete(beta_j,(0),axis=0) #intercept value in 0
        
        pip_j=np.loadtxt(JODIE_path+"mean_prob.txt") #PER SNP NOT PER VARIABLE 
        
        ### filter confidence beta: pip,min beta and var threshold
        filt_beta=beta_j
        if beta_thr:
            filt_beta=np.array(np.where(np.array(filt_beta)>=beta_thr,filt_beta,0).tolist()) #replace value with 0 if <=threshold
        if pip_thr:
            pip_j_ind=np.where(np.array(pip_j)>pip_thr)[0] #filtered pip by threshold
            
            filt_pip=np.zeros(beta_j.shape)
            for i in pip_j_ind:
                filt_pip[i,:]=filt_beta[i,:]
            filt_beta=filt_pip
            
        if var_thr:
            #variance jodie
            var=np.abs(np.diagonal(np.array(covars)))
            
            if var_comp: #the filter is applied for each component (only allow coefficients from the component that passed threshold)
                #per component
                select_comp=var>var_thr #components/parameters selected 
                filt_beta=filt_beta*select_comp
            else:
                #any variance value > threshold 0.05
                select=np.any(var>var_thr)
                if select==False: #none of the components has a significant variance
                    filt_beta=np.zeros(beta_j.shape) #all 0
            
        beta_j=np.reshape(beta_j,beta_j.shape[0]*beta_j.shape[1]) #reshape beta to plane list
        
    except:
        print("JODIE did not converged")
        return(False)
    
    return(covars,beta_j,pip_j,filt_beta) #beta=plane list of raw values, filt_beta=array shaped by parameters, filtered based on thresholds

def loadSUSIEResults(susie_path,name,beta_thr=False,pip_thr=False):
    beta=zarr.load(susie_path+"beta_"+name+".zarr")
    beta=beta[:,0]
    beta=np.delete(beta,0) #intercept 0
    
    pip=zarr.load(susie_path+"pip_"+name+".zarr") #per SNP
    pip=pip[:,0]
    
    #filter beta based on min value and pip
    filt_beta=beta
    if beta_thr:
        filt_beta=np.where(np.array(beta)>=beta_thr,filt_beta,0).tolist()
    if pip_thr:
        filt_beta=np.where(np.array(pip)>=pip_thr,filt_beta,0).tolist()
    filt_beta=np.reshape(filt_beta,[int(len(beta)/4),4])
    return(beta,pip,filt_beta)

def loadLASSOResults(path,name,beta_thr=False,pip_thr=False):
    beta=np.loadtxt(path+"meanBeta_"+name+".txt",delimiter=",")
    pip=np.loadtxt(path+"pip_"+name+".txt",delimiter=",")
    #filter beta based on min value and pip
    filt_beta=beta
    if beta_thr:
        filt_beta=np.where(np.array(beta)>=beta_thr,filt_beta,0).tolist()
    if pip_thr:
        filt_beta=np.where(np.array(pip)>=pip_thr,filt_beta,0).tolist()
    filt_beta=np.reshape(filt_beta,[int(len(beta)/4),4])
    return(beta,pip,filt_beta)

###analysis
def r2_determination(y_true, y_pred): #Ilse's code: coefficient of determination
    y_mean = np.mean(y_true)
    ssres = np.sum((y_true-y_pred)**2)
    sstot = np.sum((y_true-y_mean)**2)
    r2 = 1-ssres/sstot
    return r2

def r2_pearson(y_true,y_pred): #pearson correlation coefficient
    if (np.std(y_true)*np.std(y_pred))==0:
        r2=0
    else:
        r=np.cov(y_true,y_pred)[0,1]/(np.std(y_true)*np.std(y_pred))
        r2=r**2
    return r2

def ROC(true_beta,pred_beta):
    #true positives
    TP=true_beta & pred_beta #matrix form
    tp=np.sum(TP)
    
    #false negatives
    FN=np.logical_and(true_beta,np.logical_not(pred_beta))
    fn=np.sum(FN)
    tpr=tp/(tp+fn)
    
    #false positives
    FP=np.logical_and(pred_beta,np.logical_not(true_beta))
    fp=np.sum(FP)
    
    #true negative
    TN=np.logical_and(np.logical_not(true_beta),np.logical_not(pred_beta))
    tn=np.sum(TN)
    
    fpr=fp/(fp+tn)
    
    #false discovery rate
    if fp+tp==0:
        fdr=0
    else:
        fdr=fp/(fp+tp)
    
    return(tpr,fpr,fdr,FP,TP)

def ROC_SNP_par(true_beta,pred_beta): #function ROC per component and SNP, FD/TP individual components
    #convert boolean
    true_beta=true_beta!=0
    pred_beta=pred_beta!=0
    
    tpr,fpr,fdr,FP,TP=ROC(true_beta,pred_beta)
    
    #true positives per parameter
    tpC,tpM,tpF,tpI=np.sum(TP,axis=0)
    
    #false positives per parameter
    fpC,fpM,fpF,fpI=np.sum(FP,axis=0)
    
    ####PER SNP
    true_beta_snp=np.sum(true_beta,axis=1)!=0
    pred_beta_snp=np.sum(pred_beta,axis=1)!=0
    
    tpr_snp,fpr_snp,fdr_snp,FP_snp,TP_snp=ROC(true_beta_snp,pred_beta_snp)

    ROC_parameter=[tpr,fpr,fdr]
    ROC_snp=[tpr_snp,fpr_snp,fdr_snp]
    FP_parameter=[fpC,fpM,fpF,fpI]
    TP_parameter=[tpC,tpM,tpF,tpI]
    
    return(ROC_parameter,ROC_snp,FP_parameter,TP_parameter)

def analysisSimJODIE_SUSIE_LASSO(X,path,name,beta_thr,pip_thr,var_thr): #load results simulations per method, get r2, ROC
    #########################################################
    #load original beta
    beta_real=zarr.load(path+"out_beta_y/"+"beta_"+name+".zarr")
    y_real=zarr.load(path+"out_beta_y/"+"y_"+name+".zarr")
    
    #load JODIE
    try:
        var_j,beta_j,pip_j,filt_j=loadJODIEResults(path+"out_jodie/"+name+"/",beta_thr,pip_thr,var_thr) #var_thr
    except:
        return(False)
    
    #load SuSIE
    beta_s,pip_s,filt_s=loadSUSIEResults(path+"out_susieR/",name,beta_thr,pip_thr)
    
    #load LASSO
    beta_l,pip_l,filt_l=loadLASSOResults(path+"out_lasso/",name,beta_thr,pip_thr)
    
    ################################
    # R2
    r2_pea_j=r2_pearson(y_real, np.matmul(X,beta_j))
    r2_det_j=r2_determination(y_real, np.matmul(X,beta_j))
    
    r2_pea_s=r2_pearson(y_real, np.matmul(X,beta_s))
    r2_det_s=r2_determination(y_real, np.matmul(X,beta_s))
    
    r2_pea_l=r2_pearson(y_real, np.matmul(X,beta_l))
    r2_det_l=r2_determination(y_real, np.matmul(X,beta_l))
    
    #################################
    #ROC
    true_beta=np.reshape(beta_real,[int(len(beta_real)/4),4])
    ROC_j,ROC_snp_j,FP_j,TP_j=ROC_SNP_par(true_beta,filt_j)
    ROC_s,ROC_snp_s,FP_s,TP_s=ROC_SNP_par(true_beta,filt_s)
    ROC_l,ROC_snp_l,FP_l,TP_l=ROC_SNP_par(true_beta,filt_l)
    
    #OUTPUT
    #injected index
    #index_snp=np.where(true_beta!=0)[0].tolist()
    #maf_snp=[maf[i] for i in index_snp]
    real=[true_beta,y_real]
    #methods
    jodie=[beta_j,pip_j,var_j,filt_j,r2_pea_j,r2_det_j,ROC_j,ROC_snp_j,FP_j,TP_j]
    susie=[beta_s,pip_s,filt_s,r2_pea_s,r2_det_s,ROC_s,ROC_snp_s,FP_s,TP_s]
    lasso=[beta_l,pip_l,filt_l,r2_pea_l,r2_det_l,ROC_l,ROC_snp_l,FP_l,TP_l]
    return(real,jodie,susie,lasso)


#Paths
####################
#Test
pathx="chr22_files/x_22.zarr/"
path="sim_final/"
out_dir="sim_final/out/"

#####################
#load X
X=zarr.load(pathx)

####SET VALUES for threshold
pip_thr=0.95
beta_thr=0.01
var_thr=0.05
########################################################
### For all files
real=[];jodie=[];susie=[];lasso=[];patterns=[]
for sim in glob.glob(path+"out_beta_y/beta_"+"*.zarr"):
   name_file=sim.split("/")[len(sim.split("/"))-1].split("beta_")[1].split(".zarr")[0] #get only iteration_pattern
   pattern=name_file.split("_")[1]
 
   result=analysisSimJODIE_SUSIE_LASSO(X, path, name_file, beta_thr, pip_thr, var_thr)
   if result is False:
           print(name_file)
           noConverge.append(name_file)
   else:
       real.append(result[0]);jodie.append(result[1]);susie.append(result[2]);lasso.append(result[3]);patterns.append(pattern)
     
##############################################
#data frame
resultsSim=pd.DataFrame({"pattern":patterns,"r2_Pearson_JODIE":[x[4] for x in jodie],"r2_determination_JODIE":[x[5] for x in jodie],
                         "[TPR,FPR,FDR]parameter_JODIE":[x[6] for x in jodie],"[TPR,FPR,FDR]SNP_JODIE":[x[7] for x in jodie],"FP_JODIE":[x[8] for x in jodie],
                         "TP_JODIE":[x[9] for x in jodie],"var_JODIE":[x[2] for x in jodie],
                         "r2_Pearson_SUSIE":[x[3] for x in susie],"r2_determination_SUSIE":[x[4] for x in susie],"[TPR,FPR,FDR]parameter_SUSIE":[x[5] for x in susie],
                         "[TPR,FPR,FDR]SNP_SUSIE":[x[6] for x in susie],"FP_SUSIE":[x[7] for x in susie],"TP_SUSIE":[x[8] for x in susie],
                         "r2_Pearson_LASSO":[x[3] for x in lasso],"r2_determination_LASSO":[x[4] for x in lasso],"[TPR,FPR,FDR]parameter_LASSO":[x[5] for x in lasso],
                         "[TPR,FPR,FDR]SNP_LASSO":[x[6] for x in lasso],"FP_LASSO":[x[7] for x in lasso],"TP_LASSO":[x[8] for x in lasso]})


#save dataframe
resultsSim.to_csv(out_dir+"resultsSimulations.csv",index=False)

#save list of lists to picle file
with open(out_dir+"true.pkl", "wb") as f:
    pickle.dump(real, f)

with open(out_dir+"jodie.pkl", "wb") as f:
    pickle.dump(jodie, f)
    
with open(out_dir+"susie.pkl", "wb") as f:
    pickle.dump(susie,f)

with open(out_dir+"lasso.pkl", "wb") as f:
    pickle.dump(lasso, f)

resultsSim.to_pickle(out_dir+"resultsSim.pkl")
