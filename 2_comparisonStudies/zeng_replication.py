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
import seaborn as sns
import math

### jodie
path_jodie="POE_mQTLs/outputs/poe_probes.csv"
jodie=pd.read_csv(path_jodie)

### gtca
path="h2_siblings/2zeng/2run/out/"

#sibs variance is g,k,s
var_sibs=[]
chro=[];probes=[];no_converged=[]
for probe in next(os.walk(path))[1]:
    path_file=path+probe+"/sibs_const.hsq" #/sp_const.hsq
    try:
        df=pd.read_csv(path_file,sep="\t")  # or sep="\s+" if space-separated
        sibs_var=df.loc[df["Source"] == "V(G3)", "Variance"].iloc[0]
        
        var_sibs.append(sibs_var)
    
        chro.append(int(probe.split("_")[0]))
        probes.append(int(probe.split("_")[1]))
    except:
        no_converged.append(probe)

gcta_df=pd.DataFrame({"Chromosome":chro,"Index":probes,"varSib":var_sibs})

### merge df
compare_df=gcta_df.merge(jodie[["Chromosome","Index","Probe ID","varI"]],on="Index")

######################################################################
def loadBestFit(path,models,const):
    LRT=0
    for model in models:
        path_file=path+"/"+model+const+".hsq"
        try:
            df=pd.read_csv(path_file,sep="\t")
            LRT_model=df.loc[df["Source"]=="LRT","Variance"].iloc[0]
          
        except: #not converged
            LRT_model=0
        
        if LRT_model>LRT:
              LRT=LRT_model
              best_fit=model
              pval=df.loc[df["Source"]=="Pval","Variance"].iloc[0]  

    if LRT==0: #not converged for all models
        best_fit=np.nan
        pval="no converged"
        
    return(LRT,best_fit,pval)

def loadVariances(path,best_fit,const):
    path_file=path+"/"+best_fit+const+".hsq"#path+probe+"/"+best_fit+"_const.hsq"
    df=pd.read_csv(path_file,sep="\t")
    
    varG=df.loc[df["Source"]=="V(G1)","Variance"].iloc[0]
    varK=df.loc[df["Source"]=="V(G2)","Variance"].iloc[0]
    if best_fit=="null":
        varS=0;varSm=0;varSp=0
    elif best_fit=="sibs":
        varS=df.loc[df["Source"]=="V(G3)","Variance"].iloc[0]
        varSm=0;varSp=0;
    elif best_fit=="sm":
        varS=0;varSp=0
        varSm=(df.loc[df["Source"]=="V(G3)","Variance"].iloc[0])
    elif best_fit=="sp":
        varS=0;varSm=0
        varSp=(df.loc[df["Source"]=="V(G3)","Variance"].iloc[0])
    elif best_fit=="smp":
        varS=df.loc[df["Source"]=="V(G3)","Variance"].iloc[0]
        varSm=0;varSp=0
    return(varG,varK,varS,varSm,varSp)

#### Getting LRT from all probes, getting variances from best fit and add column of best fit model
## chr probe model varG varK varS varSm varSp 
p_threshold=0.05

# best fit with LRT in both constrained and unconstrained with Pvalue of the selected being less than threshold
models=["null","sibs","sm","sp","smp"]#

fit_no=[];fit_const=[];fit_models=[];rep=[]
pvals=[];pvals_no=[];LRT=[];LRT_no=[]
chroms=[];probes=[];varG=[];varK=[];varS=[];varSm=[];varSp=[]
status=[];

for probe in next(os.walk(path))[1]:
    #constrained
    LRT_const,best_fit,pval=loadBestFit(path+probe, models, "_const")
    
    #not constrained
    LRT_no_c,best_fit_no,pval_no=loadBestFit(path+probe, models, "_no")
    
    fit_no.append(best_fit_no);fit_const.append(best_fit)
    LRT.append(LRT_const);LRT_no.append(LRT_no_c)
    pvals.append(pval);pvals_no.append(pval_no)
    
    ####
    if best_fit==best_fit_no: #constrained and unconstrained  match
        fit_models.append(best_fit)
        status.append("consistent")
        #if best_fit=="null"status.append("null")
        varG_fit,varK_fit,varS_fit,varSm_fit,varSp_fit=loadVariances(path+probe,best_fit,"_const")
        
    else: #no matching constrained and unconstrained
        
        if best_fit!="null" and (best_fit_no=="null" or best_fit_no=="no converged") and LRT_const>LRT_no_c:
            fit_models.append(best_fit)
            status.append("const. improvement")
            varG_fit,varK_fit,varS_fit,varSm_fit,varSp_fit=loadVariances(path+probe,best_fit,"_const")
            
        elif (best_fit=="null" or best_fit=="no converged") and LRT_no_c>=LRT_const: #unconstrained better LRT compared to null 
            fit_models.append(best_fit_no) 
            status.append("recovered")
            varG_fit,varK_fit,varS_fit,varSm_fit,varSp_fit=loadVariances(path+probe,best_fit_no,"_no")
            
        elif best_fit!="null": #take constrained if unconstrained = null
            fit_models.append(best_fit)
            status.append("costrained")
            varG_fit,varK_fit,varS_fit,varSm_fit,varSp_fit=loadVariances(path+probe,best_fit,"_const")
            
        else:
            print(best_fit,LRT_const,best_fit_no,LRT_no_c)
       
    chroms.append(int(probe.split("_")[0]))
    probes.append(int(probe.split("_")[1]))
    varG.append(varG_fit);varK.append(varK_fit);varS.append(varS_fit);varSm.append(varSm_fit);varSp.append(varSp_fit)
        
        
df_fits=pd.DataFrame({"Chromosome":chroms,"Index":probes,"model":fit_models,"status":status,"varG":varG,
                      "varK":varK,"varS":varS,"varSm":varSm,"varSp":varSp,"pval":pvals,"LRT":LRT,
                      "model constrained":fit_const,"model no":fit_no,"LRT no":LRT_no,
                      "pval no":pvals_no})


df_fits=df_fits.merge(jodie[["Chromosome","Index","Probe ID","varC","varM","varF","varI","covCI","covCM","covCF","covMF","covMI","covFI","Gene"]],
                      on=["Chromosome","Index"])



plt.figure()
plt.hist(df_fits["model"])

plt.figure()
plt.hist(df_fits["status"])


a=df_fits[["model","varI","varC","varM","varF","varG","varK","varS","varSm","varSp"]]

fits=df_fits.copy()
fits=fits.dropna(subset="model") #only replicated

fits.to_csv("h2_siblings/2zeng/3analysis/fits.csv",index=False)
df_fits.to_csv("h2_siblings/2zeng/3analysis/zeng_rep.csv",index=False)
