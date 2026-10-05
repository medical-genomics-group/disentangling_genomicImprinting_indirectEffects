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

###################################################################################
# overlap previous studies (Rosenski, Cuellar and Zeng)
##################################################################################
buffer=2000000
path_save="POE_mQTLs/outputs/"

#### LOAD DATA
################################### # Rosenski 
# Data S2. Known iDMRs with ubiquituous hyper/hypo methylated sub-regions removed.
known_impt_dir="data/knownImp_coordinates.xlsx"
known_imp=pd.read_excel(known_impt_dir,skiprows=[0,1],sheet_name=0)
known_imp["Chr"]=[int(x.split("chr")[1]) for x in known_imp["chrom"].tolist()]
known_imp["Gene"]=[x.split(":")[0] for x in known_imp["ICR name"]]
known_imp=known_imp[["Chr","start","end","Gene"]]
known_imp["Source"]="Rosenski iDMR"
known_imp["Position"]=(known_imp["end"]-known_imp["start"])/2+known_imp["start"] ###middle

####Data S4. Parental-ASM SNP loci and covering bimodal region. Shown are p-values (Fisher's exact), as well as Benjamini-Hochberg adjusted FDR values.
path_ros="data/S4_SNPs_parental_ASMs.xlsx"
ros=pd.read_excel(path_ros,skiprows=2)
ros=ros[["chrom",'SNP pos (hg19)','dbSNP id (gnomAD 2.1.1 asj)','ICR known associated genes',"Parental ASM region (hg19)"]]
ros["Chr"]=[int(x.split("chr")[1].split(":")[0]) for x in ros["Parental ASM region (hg19)"].tolist()]
ros["start"]=[int(x.split(":")[1].split("-")[0]) for x in ros["Parental ASM region (hg19)"].tolist()]
ros["end"]=[int(x.split(":")[1].split("-")[1]) for x in ros["Parental ASM region (hg19)"].tolist()]

ros=ros.rename(columns={"SNP pos (hg19)":"Position SNP","dbSNP id (gnomAD 2.1.1 asj)":"RSID",
                        "ICR known associated genes":"Gene"})
ros=ros[["Chr","start","end","RSID","Position SNP","Gene"]]
ros["Source"]="Rosenski p.ASM"
ros["Position"]=(ros["end"]-ros["start"])/2+ros["start"] ###middle

####################################################
###Cuellar Partida
path_parS1="data/Partida_S1.xlsx"
parS1=pd.read_excel(path_parS1) 
#Supplementary Table 1.Supplementary Table 1. Genome-wide statistically significant parent-of-origin effects on DNA methylation. Only results for the top associated SNP at each CpG site are included. Results where the test of association did not reach nominal significance (p-value >0.05) were not stored.																		
columns=list(parS1.columns.values)
keep_columns=columns[0:3]+["SNP"]
parS1=parS1[["Chr","CpG_BP","CpG_ID","Gene","SNP"]] #probes and genes
parS1=parS1.rename(columns={"CpG_BP":"Position","CpG_ID":"Probe ID","SNP":"RSID"})
parS1["Source"]="Cuellar POE-CpG"

###known imprinted Loci Cuellar Partida S2
#Supplementary Table 2. List of loci implicated in previous studies of imprinting. A locus was defined to be at least 2Mb apart from one another. Transcription start and end positions are based on the GRCh37 Human genome assembly.																				
path_parS2="data/Partida_S2_knownImprinted.xlsx"
parS2=pd.read_excel(path_parS2)
parS2=parS2[["Chr","Start","End","Gene"]]
parS2=parS2.rename(columns={"Start":"start","End":"end"})
parS2["Position"]=(parS2["end"]-parS2["start"])/2+parS2["start"] ##middle
parS2["Source"]="Cuellar Imprint"

#######################################
###### Zeng
zeng_poe="data/Zeng_POE_cpg.xlsx"
zeng_poe=pd.read_excel(zeng_poe)
zeng_poe=zeng_poe.rename(columns={"CpG":"Probe ID","CHR":"Chr"})
zeng_poe["Source"]="Zeng"
######################################

### load my data
my_probes=pd.read_csv("POE_mQTLs/outputs/poe_probes.csv") 
my_probes=my_probes[["Chromosome","Probe ID","Position","Gene"]]
my_probes=my_probes.rename(columns={"Chromosome":"Chr"})
my_probes["Source"]="Ours"

def overlapStudies(list_studies):
    probes_merged=pd.concat(list_studies,ignore_index=True) 

    ###how many by location?
    #position within buffer +/- 
    for chro in np.unique(probes_merged["Chr"]).tolist():
        probes_merged_chr=probes_merged[(probes_merged["Chr"]==chro)] #subset by chr and not overlap already by id
        # & (probes_merged["Overlap"]!="ID")
        my_probes_chr=probes_merged_chr[probes_merged_chr["Source"]=="Ours"]
        compare=probes_merged_chr[probes_merged_chr["Source"]!="Ours"]
        
        pos=my_probes_chr.Position.values
        end=np.array([x+buffer for x in compare["Position"]])
        start=np.array([x-buffer for x in compare["Position"]])
        
        i,j=np.where((pos[:,None]>=start)&(pos[:,None]<=end))
        i=np.unique(i).tolist();j=np.unique(j).tolist()
        if len(i)>0:
            # Map back to original dataframe indices
            idx_my  = my_probes_chr.iloc[i].index #index of chr lost actual index in p_mer_chr
            idx_cmp = compare.iloc[j].index

            probes_merged.loc[idx_my,"Overlap"]="Position"
            probes_merged.loc[idx_cmp,"Overlap"]="Position"
        
    id_dup=probes_merged.duplicated(subset="Probe ID",keep=False).tolist()
    not_nan=probes_merged["Probe ID"].isna().tolist()
    not_nan=np.bitwise_not(not_nan).tolist() #negate
    probes_merged["Overlap"]=["ID" if (x and y) else i for x,y,i in zip(id_dup,not_nan,probes_merged["Overlap"].tolist())]#id duplicated but not nan

    ###how many by id?
    same_probe=probes_merged[(probes_merged["Source"]=="Ours")&(probes_merged["Overlap"]=="ID")]
    print("Number of probes replicated by probe id: ",len(same_probe))

    ###how many by distance
    print("Number of probes replicated by distance <=",buffer,": ",
          len(probes_merged[(probes_merged["Source"]=="Ours")&(probes_merged["Overlap"]=="Position")])) 

    ##double check replicated
    replicated=probes_merged[(probes_merged["Source"]=="Ours")&(probes_merged["Overlap"].notnull())]

    novel=probes_merged[(probes_merged["Source"]=="Ours")&(probes_merged["Overlap"].isna())] 
    print("Number of novel POE dependent probes: ",len(novel)) #42

    return(probes_merged,replicated,novel)

print("Including zeng: ")
merged_z,replicated_z,novel_z=overlapStudies([my_probes,known_imp,ros,parS1,parS2,zeng_poe]) #including zeng
print("Only external cohorts")
merged,replicated,novel=overlapStudies([my_probes,known_imp,ros,parS1,parS2]) #excluding zeng


##save dfs
merged.to_csv(path_save+"comparisonStudies.csv",index=False)


########################################
#all probes screen 2 
all_probes=pd.read_csv("POE_mQTLs/outputs/screen2.csv")

#1. novel variances
novel_var=all_probes[all_probes["Probe ID"].isin(novel["Probe ID"])]
print("Nomber of novel probes: ",len(novel_var))

#2. replicated
replicated=all_probes[all_probes["Probe ID"].isin(replicated["Probe ID"])]
print("Number of probes in the replicated category: ",len(replicated)) 

#3. other studies by id 
otherStudies=merged[(merged["Source"]!="Ours")&(merged["Overlap"].isna())]
otherStudies=otherStudies.dropna(subset=["Probe ID"]) 

#from which we have data for
otherStudies=all_probes[all_probes["Probe ID"].isin(otherStudies["Probe ID"].tolist())] 
print("Number of probes in other studies category: ",len(otherStudies))


novel.to_csv(path_save+"novel_studiesComparison.csv")
replicated.to_csv(path_save+"overlap_studiesComparison.csv")
otherStudies.to_csv(path_save+"others_studiesComparison.csv")

####
# Genes
###
parS1_genes=parS1.Gene.dropna().values.tolist()
parS2_genes=parS2.Gene.dropna().values.tolist()
par_genes=np.unique(parS1_genes+parS2_genes).tolist()

known_genes=known_imp.Gene.dropna().values.tolist()
ros_genes_raw=ros.Gene.dropna().values.tolist()
ros_genes=[]
replacement=({"Multiple":"","transcripts":""," and microRNA cluster":""," ":""})
for row in ros_genes_raw:
    for old,new in replacement.items():
        row=row.replace(old,new)
    row=row.split(",")
    ros_genes.append(row)
    
ros_genes=[x for xs in ros_genes for x in xs]
ros_genes=np.unique(ros_genes).tolist()

studies_genes=par_genes+known_genes+ros_genes
studies_genes=np.unique(studies_genes).tolist()

with open(path_save+'genes_previousStudies.txt', 'w') as f:
    for line in studies_genes:
        f.write(f"{line}\n")


#load all genes ours
ourGenes_path="POE_mQTLs/outputs/OurGenes.txt"
ourGenes=np.loadtxt(ourGenes_path,dtype=str).tolist()
overlap=list(set(studies_genes).intersection(set(ourGenes))) 
novel_genes=list(set(ourGenes)-set(studies_genes)) 

with open(path_save+'genes_overlap.txt', 'w') as f:
    for line in overlap:
        f.write(f"{line}\n")

with open(path_save+'genes_novel.txt', 'w') as f:
    for line in novel_genes:
        f.write(f"{line}\n")


