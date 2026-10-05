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

################################################################################
## downstream associations eQTLs, pQTLs, QTLs others
clusters_df_path="LD_fineMap/4run_regression/mQTL_clusterAnnotate.csv"
my_snps=pd.read_csv(clusters_df_path)

path_eqtls="eQTLs/out/"
path_save="POE_mQTLs/outputs/"

ld_df="LD_fineMap/5matchSNPs/r2_mQTL_fineMap.csv"
ld_df=pd.read_csv(ld_df)
ld_df=ld_df.rename(columns={"SNP_B":"RSID","SNP_A":"Fine-mapped"})
ld_df=pd.merge(ld_df,my_snps[["Probe ID","RSID"]],left_on="Fine-mapped",right_on="RSID")
ld_df=ld_df.rename(columns={"RSID_x":"RSID"})

pval_thr=0.005
################ ONEK1K
onek1k_path="eQTLs/esnp_table.tsv"
onek1k=pd.read_csv(onek1k_path,sep="\t")
onek1k=onek1k[onek1k["P_VALUE"]<=pval_thr]

## hits onek1k
mqtls_onek1k=pd.merge(my_snps,onek1k,on="RSID") #7
mqtls_onek1k.to_csv(path_eqtls+"matchID_onek1k",index=False) 
print("Matches onek1k cluster 1Mb selected snps: ",len(mqtls_onek1k))


onek1k_matches=pd.merge(ld_df[["RSID","Fine-mapped","R2","Probe ID"]],onek1k,on="RSID") 
print("Matches with onek1k done!")
print("Number of matches with clump: ",len(onek1k_matches)) #185
onek1k_matches.to_csv(path_save+"matches_clump_onek1k.csv",index=False)

onek1k_snps=onek1k_matches.drop(columns=["Fine-mapped","R2","Probe ID"])
onek1k_snps=onek1k_snps.drop_duplicates()
print("number of unique SNPs matches: ",len(onek1k_snps),"\n")

onek1k_probes=set(onek1k_matches["Probe ID"].tolist())
print("number of unique probes matches: ",len(onek1k_probes),"\n")
print("/////////////////////////\n") 

##unique genes matched
genes_onek1k=onek1k_matches.iloc[:,6].tolist()
genes_onek1k=np.unique(genes_onek1k).tolist()

with open(path_save+'onek1k_genes.txt','w') as f:
    for line in genes_onek1k:
        f.write(f"{line}\n")

#########################
### GTEx eQTLs
#######################
gtex_path="eQTLs/GTEx/Whole_Blood.v11.eGenes.txt.gz"
gtex_eqtls=pd.read_csv(gtex_path,delimiter="\t")
#gtex_eqtls=gtex_eqtls[gtex_eqtls["qval"]<=0.05] #list of eGenes
gtex_eqtls=gtex_eqtls[gtex_eqtls["pval_nominal"]<=pval_thr]
gtex_eqtls=gtex_eqtls.rename(columns={"rs_id_dbSNP157_GRCh38p14":"RSID"})

#overlap
mqtls_gtex=pd.merge(my_snps,gtex_eqtls,on="RSID") #0 no overlap by rsid directly
print("Matches gtex cluster 1Mb selected snps: ",len(mqtls_gtex))
#by clump
gtex_matches=pd.merge(ld_df[["RSID","Fine-mapped","R2","Probe ID"]],gtex_eqtls,on="RSID") 
print("Matches with gtex done!")
print("Number of matches with clump: ",len(gtex_matches)) #185
gtex_matches.to_csv(path_save+"matches_clump_gtex.csv",index=False)

gtex_snps=gtex_matches.drop(columns=["Fine-mapped","R2","Probe ID"])
gtex_snps=gtex_snps.drop_duplicates()
print("number of unique SNPs matches: ",len(gtex_snps),"\n")

gtex_probes=set(gtex_matches["Probe ID"].tolist())
print("number of unique probes matches: ",len(gtex_probes),"\n")
print("/////////////////////////\n") 

####################
#PQTLS
##################
#Hofmeister
qtls_hof=pd.read_excel("eQTLs/hofmeister_POEs.xlsx") 
QTLsHof_mQTLs=pd.merge(qtls_hof,my_snps,left_on="SNP ID",right_on="RSID") 
print("Matches hofmeister qtls cluster 1Mb selected snps: ",len(QTLsHof_mQTLs))
#by clump
hof_qtls_matches=pd.merge(ld_df[["RSID","Fine-mapped","R2","Probe ID"]],qtls_hof,right_on="SNP ID",left_on="RSID")
print("Matches with Hofmeister qtls done!")
print("Number of matches: ",len(hof_qtls_matches)) 
hof_qtls_matches.to_csv(path_save+"matches_clump_hofQTLs.csv",index=False)
print("CSV file saved")

hof_qtls_snps=hof_qtls_matches.drop(columns=["Fine-mapped","R2","Probe ID"])
hof_qtls_snps=hof_qtls_snps.drop_duplicates()
print("number of unique SNPs matches: ",len(hof_qtls_snps),"\n")

hof_qtls_probes=set(hof_qtls_matches["Probe ID"].tolist())
print("number of unique probes matches: ",len(hof_qtls_probes),"\n")
print("/////////////////////////\n") 


pQTLs_hof=pd.read_excel("eQTLs/pQTLs_Hofmeister.xlsx") 
pQTLsHof_mQTLs=pd.merge(pQTLs_hof,my_snps,on="RSID") 
print("Matches hofmeister pQTLs cluster 1Mb selected snps: ",len(pQTLsHof_mQTLs))

#by clump
hof_pQTLs_matches=pd.merge(ld_df[["RSID","Fine-mapped","R2","Probe ID"]],pQTLs_hof,on="RSID")
print("Matches with Hofmeister pQTLs done!") 
print("Number of matches: ",len(hof_pQTLs_matches))
hof_pQTLs_matches.to_csv(path_save+"matches_clump_hofpQTLs.csv",index=False)
print("CSV file saved")

hof_pQTLs_snps=hof_pQTLs_matches.drop(columns=["Fine-mapped","R2","Probe ID"])
hof_pQTLs_snps=hof_pQTLs_snps.drop_duplicates()
print("number of unique SNPs matches: ",len(hof_pQTLs_snps),"\n")

hof_pQTLs_probes=set(hof_pQTLs_matches["Probe ID"].tolist())
print("number of unique probes matches: ",len(hof_pQTLs_probes),"\n")
print("/////////////////////////\n") 


###UK biobank
###
path_uk="pQTLs/UK_pQTLs_summary.xlsx"
uk_pqtls=pd.read_excel(path_uk,sheet_name="ST9",skiprows=[0,1,2,3])
uk_pqtls=uk_pqtls[uk_pqtls["log10(p) (discovery)"]>2.301] #p<0.005

pQTLs_UK_mQTLs=pd.merge(uk_pqtls,my_snps,left_on="rsID",right_on="RSID") 
print("Matches UK pQTLs cluster 1Mb selected snps: ",len(pQTLs_UK_mQTLs))

#by clump
uk_pQTLs_matches=pd.merge(ld_df[["RSID","Fine-mapped","R2","Probe ID"]],uk_pqtls,left_on="RSID",right_on="rsID")
print("Matches with UK pQTLs done!") 
print("Number of matches: ",len(uk_pQTLs_matches))
uk_pQTLs_matches.to_csv(path_save+"matches_clump_UKpQTLs.csv",index=False)
print("CSV file saved")

uk_pQTLs_snps=uk_pQTLs_matches.drop(columns=["Fine-mapped","R2","Probe ID"])
uk_pQTLs_snps=uk_pQTLs_snps.drop_duplicates()
print("number of unique SNPs matches: ",len(uk_pQTLs_snps),"\n")

uk_pQTLs_probes=set(uk_pQTLs_matches["Probe ID"].tolist())
print("number of unique probes matches: ",len(uk_pQTLs_probes),"\n")
print("/////////////////////////\n") 

