"""
GENES
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
import pickle


path_out="POE_mQTLs/outputs/"

clusters_df_path="LD_fineMap/4run_regression/mQTL_clusterAnnotate.csv"
my_snps=pd.read_csv(clusters_df_path)

my_snps["RSID"]=my_snps["rsID cluster SNP"]
my_snps=my_snps.drop("rsID cluster SNP",axis=1)
my_snps=my_snps.drop_duplicates(["RSID"])

my_probes_path="POE_mQTLs/outputs/poe_probes.csv"
my_probes=pd.read_csv(my_probes_path)

##### create list of genes, ensembl + probes
#https://gdc.cancer.gov/content/improved-dna-methylation-array-probe-annotation
#https://zwdzwd.github.io/InfiniumAnnotation#reference
#Human MSA - Goldberg et al. Scalable Screening of Ternary-Code DNA Methylation Dynamics Associated with Human Traits, Cell Genomics 2025
msa_path="genes/MSA.hg38.manifest.gencode.v41.tsv.gz"
msa=pd.read_csv(msa_path,delimiter="\t")
msa=msa.dropna(subset=["CpG_chrm","genesUniq"])
msa["Chromosome"]=[x.split("chr")[1] for x in msa["CpG_chrm"].tolist()]
msa["Probe ID"]=[x.split("_")[0] for x in msa["probeID"].tolist()]
msa=msa[(msa["Chromosome"]!="Y") & (msa["Chromosome"]!="M")]
msa=msa.replace("X","23")
msa["Chromosome"]=[int(x) for x in msa["Chromosome"].tolist()]

#merged with my probes
msa_poe=pd.merge(msa[["Chromosome","Probe ID","genesUniq","transcriptTypes"]],my_probes,
                 on=["Chromosome","Probe ID"])

msa_poe.to_csv("POE_mQTLs/outputs/msaPOE_genes.csv",index=False)

ensembl_path="POE_mQTLs/outputs/mQTLs_ensembl.txt"
ensembl=pd.read_csv(ensembl_path,delimiter="\t")
ensembl=ensembl.dropna(subset="SYMBOL")

ensembl_genes=ensembl[["#Uploaded_variation","SYMBOL"]]
ensembl_genes=ensembl_genes[ensembl_genes["SYMBOL"]!="-"]

ensembl_genes=ensembl_genes.drop_duplicates()

snps_genes=pd.merge(my_snps[["Chromosome","Index probe","Probe ID","RSID","Gene Probe"]], #"Imprinting pattern","Non imprinting effects"
                    ensembl_genes,left_on="RSID",right_on="#Uploaded_variation",how="left")

def col2list(df,x_cols):
    agg_dict = {
        col: lambda x: ";".join(dict.fromkeys(map(str, x)))
        for col in x_cols
    }

    df_out=df.groupby(
        [c for c in df.columns if c not in x_cols],
        as_index=False,dropna=False
    ).agg(agg_dict)
    return(df_out)

snps_genes=col2list(snps_genes,["SYMBOL"])

snps_genes=pd.merge(snps_genes,msa[["Chromosome","Probe ID","genesUniq"]],on=["Chromosome","Probe ID"],how="left")
snps_genes=snps_genes.drop_duplicates()
snps_genes=snps_genes.rename(columns={"Gene Probe":"GenScot","SYMBOL":"ENSEMBL","genesUniq":"MSA"})
snps_genes.drop("#Uploaded_variation",axis=1,inplace=True)

snps_genes.to_csv("POE_mQTLs/outputs/SNPs_genes.csv",index=False)


probes_genes_dbs=pd.merge(my_probes[["Chromosome","Index","Probe ID","Gene"]],
                          msa[["Chromosome","Probe ID","genesUniq","transcriptTypes"]],on=["Chromosome","Probe ID"],how="left")

probes_genes_dbs=probes_genes_dbs.drop_duplicates()

probes_genes_dbs=pd.merge(probes_genes_dbs,my_snps[["Chromosome","Probe ID","RSID"]],on=["Chromosome","Probe ID"],how="left")
probes_genes_dbs=pd.merge(probes_genes_dbs,ensembl_genes,left_on="RSID",right_on="#Uploaded_variation",how="left")
probes_genes_dbs.drop("#Uploaded_variation",axis=1,inplace=True)

probes_genes_dbs=col2list(probes_genes_dbs,["RSID","SYMBOL"])
probes_genes_dbs=probes_genes_dbs.rename(columns={"Gene":"GenScot","SYMBOL":"ENSEMBL","genesUniq":"MSA"})

######
list_snps_genes=[]
for i,row in snps_genes.iterrows():
    gene_snp=row["ENSEMBL"]
    if isinstance(row["MSA"],str):
        list_msa=row["MSA"].split(";")
    else:
        list_msa=[""]
    if isinstance(row["GenScot"],str):
        list_estonia=row["GenScot"].split(";")
    else:
        list_estonia=[""]
    if gene_snp in list_msa or gene_snp in list_estonia:
        list_snps_genes.append(gene_snp)
        
list_snps_genes=np.unique(list_snps_genes).tolist() #list of genes from ensembl seen also in annotation probes dataset

#list msa
list_msa=[]
for row in msa_poe["genesUniq"].tolist():
    for i in row.split(";"):
        list_msa.append(i)

list_msa=np.unique(list_msa).tolist()

#list annotation  genscot
list_scot=[]
for row in my_probes["Gene"].dropna().tolist():
    for i in row.split(";"):
        list_scot.append(i)
    
list_scot=np.unique(list_scot).tolist()

###merged list
merged_genes_list=list_snps_genes+list_msa+list_scot
merged_genes_list=np.unique(merged_genes_list).tolist()

with open(path_out+'allGenes.txt', 'w') as f:
    for line in merged_genes_list:
        f.write(f"{line}\n")

s1,s2,s3=set(list_msa),set(list_scot),set(list_snps_genes)
all_names=sorted(s1|s2|s3)
genes_source=pd.DataFrame({"Gene":all_names,"GenScot":[x in s2 for x in all_names],
                           "MSA":[x in s1 for x in all_names],"ENSEMBL":
                               [x in s3 for x in all_names]})
    
genes_source.to_csv(path_out+"genes_source.csv",index=False)

np.savetxt(path_out+"OurGenes.txt",genes_source.Gene.values,fmt="%s")

### ensembl
genes_probes= (
    my_probes["Gene"]
    .dropna()                
    .astype(str)              
    .str.split(";")           
    .explode()               
    .str.strip()             
    .tolist()
)

genes_probes=np.unique(genes_probes).tolist() 

#genes from ensembl rsids of mqtls 
ensembl=ensembl.replace("-",np.nan)
genes_en=(ensembl["SYMBOL"].dropna().tolist())
genes_en=np.unique(genes_en).tolist() 

#joining
genes_all=genes_en+genes_probes
genes_all=np.unique(genes_all).tolist() 

with open(path_out+'allGenes.txt', 'w') as f:
    for line in genes_all:
        f.write(f"{line}\n")
