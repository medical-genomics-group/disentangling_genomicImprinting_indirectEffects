"""
gets dataframe from all information of first screening CI model

1.get variance for all probes
2.normalized variance
3.save df
4.save information probes that weren't calculated
"""

import numpy as np
import ast
import matplotlib.pyplot as plt
import pandas as pd
import pyreadr
import zarr

normVal=3.421781583211533
var_thr=0.05

def getVariance(log):
    ind=log.index("mean_V=array")
    
    mean_V_log=log[ind:ind+74]
    beg=mean_V_log.index("[[")
    end=mean_V_log.index("]]")
    numbers=mean_V_log[beg:end+2]
    
    array_data = ast.literal_eval(numbers)
    variances= np.array(array_data)
    
    varC=variances[0,0]
    varI=variances[1,1]
    covCI=variances[0,1]
    
    varC_std=varC/normVal
    varI_std=varI/normVal
    
    ###error
    ind_error=log.index("var_V=array")
    error_V_log=log[ind_error:ind_error+120]
   
    beg=error_V_log.index("[[")
    end=error_V_log.index("]]")
    numbers=error_V_log[beg:end+2]
    
    var_V=ast.literal_eval(numbers)
    var_V=np.array(var_V)
    var_vC=var_V[0,0]
    var_vI=var_V[1,1]
    var_vCI=var_V[0,1]
    
    return(varC,varI,covCI,var_vC,var_vI,var_vCI,varC_std,varI_std)

def getLocation(chro,probe,path):
    path_name_probes=path+chro+".csv"
    probes=np.genfromtxt(path_name_probes,delimiter=",",dtype=object)
    
    name=probes[probe,0].decode()
    loc=probes[probe,1].decode()
    
    return(name,loc)    

path_name_probes="screen_full/X_matrix/out_methylation/probes_chr_"
path_chr="screenCI/values_chr.txt"
path_logs="screenCI/logs/"
path_save="outputs/"
path_name_probes="screen_full/X_matrix/out_methylation/probes_chr_"
annot_df=pyreadr.read_r("GSM-preprocessed/EPIC_AnnotationObject_df.rds")[None]
annot_df=annot_df.rename(columns={"Name":"Probe ID"})


chr_values=np.loadtxt(path_chr)
probes=chr_values[:,2]

chro=list(range(1,24))

passed=np.zeros([len(chro),2])
passed[:,0]=chro

filt=np.zeros([len(chro),3])
filt[:,0]=chro

errors=np.zeros([len(chro)])

varC=[];varI=[];covCI=[];index=[];chromosomes=[]
varC_raw=[];varI_raw=[]
var_vC=[];var_vI=[];var_vCI=[]
pos=[];names=[]


df=pd.DataFrame()

for i in range(0,len(chro)):#len(chro)
    chrom=str(chro[i])
    print("Working in chr ",chrom)
    error=[]
    missing=[]
    
    for probe in range(0,int(probes[i])): #int(probes[i]) int(probes[i]) 
        try: #try to read the file
            log_file=open(path_logs+"probe_"+chrom+"_"+str(probe)+".log",'r')
            log=log_file.read()
        except: #non existing file
            missing.append(probe)
            errors[i]+=1
            continue
        
        if ("mean_V=array" in log): #jodie gave out values  return(varC,varI,covCI,var_vC,var_vI,var_vCI,varC_std,varI_std)
            C,I,CI,eC,eI,eCI,C_std,I_std=getVariance(log)
            
            varC.append(C_std);varI.append(I_std);covCI.append(CI)
            var_vC.append(eC);var_vI.append(eI);var_vCI.append(eCI)
            varC_raw.append(C),varI_raw.append(I)
            #get location of probe
            name,loc=getLocation(chrom,probe,path_name_probes)
            index.append(probe);chromosomes.append(chrom);names.append(name);pos.append(loc)

            if I>=var_thr and np.all(np.array([var_vC,var_vI,var_vCI])<normVal): #filter significance I component
                passed[i,1]+=1
            else:
                filt[i,2]+=1 #filter by variance
              
        else:
            if ("STOPPED BY 0s IN VARIANCE AFTER 100 ITERATIONS" in log): #filt 0s
                filt[i,1]+=1 #filter by 0s in 100 iterations
            else:
                errors[i]+=1
                error.append(probe)
 
    ###
    print("Probes with errors: ",error)
    print("Probes missing: ",missing)
    

##made sure all numbers make sense
suma=errors+passed[:,1]+filt[:,1]+filt[:,2]
print("EQUAL  PROBES  SUM  ERRORS  PASSED  FILTERED 0s  Filtered Var's ")
for i in range(0,len(suma)):
    if suma[i] == probes[i]:
        print(int(passed[i,0]),"yes")
    else:
        print(int(passed[i,0]),"No:",int(probes[i]),int(suma[i]),int(errors[i]),int(passed[i,1]),int(filt[i,1]),int(filt[i,1]))
            

### save dataframe
df=pd.DataFrame({"Chromosome":chromosomes,"Index":index,"Probe ID":names,"Position":pos,"varC":varC,"varI":varI,"covCI":covCI,
                 "var varC":var_vC,"var varI":var_vI,"var covCI":var_vCI,"varC raw":varC_raw,"varI raw":varI_raw})
##annotate
df=pd.merge(df,annot_df[["UCSC_RefGene_Name","Probe ID"]],on="Probe ID",how="inner")
df=df.rename(columns={'UCSC_RefGene_Name':'Gene'})

df.to_csv(path_save+"screen1CI.csv",index=False)
print("Saved summary data frame")

### passed to second filter
second_screen=df.loc[(df["varI"]>var_thr) & (df['var varC'] <normVal) & (df['var varI'] <normVal)
                     & (df['var covCI'] <normVal)]

second_screen.to_csv(path_save+"probes_pased1.csv",index=False)

print("Probes for second screening saved")






