import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import zarr
import os
import sys
import csv

out_dir="screen_full/X_matrix/out_methylation/"

path_met="methylationData/std_methylation_chr"
path_map="methylationData/map/annotations.csv"
path_probes="GSM-preprocessed/probes_id_zarr/"

ids=pd.read_csv("methylationData/map/covars.csv") #individuals ids from probes
idLink=pd.read_csv("1.matrix/Y/idlink.csv")

trios=np.loadtxt("screen_full/trios_cMeth.ped")
childs=trios[:,0]

map_all=np.loadtxt(path_map,delimiter=",", dtype=object,usecols=(0,1,2,3,4,5,6,7)) #location methylation probes all chr

#prepare methylation data
map_all[:,3]=[x.replace('"','') for x in map_all[:,3]] #eliminate " " from string name
map_all[:,0]=[x.replace('"','') for x in map_all[:,0]] 
map_all[1:,0]=[x.split("r")[1] for x in map_all[1:,0]] #eliminate chr on rest of chr entries

#filter methylation if id is child in trio
ids=ids.rename(columns={"Sample_Sentrix_ID": "X"})
ids=ids.merge(idLink,on='X')

for chrom in range(1,23):
    chrom=str(chrom)
    map_chr=map_all[np.where(map_all[:,0]==chrom)]
    met=zarr.load(path_met+chrom+".zarr") #raw file all individuals all probes chrom
    probes=np.loadtxt(path_probes+"ids_chr"+chrom+".txt",delimiter=",", dtype=str) #probes ids

    met_child=[]
    for child in childs:
        ind=np.where(ids==child)[0][0]
        met_child.append(met[ind,:])
    
    met_child=np.array(met_child)
    zarr.save(out_dir+"methylation_"+chrom+'.zarr',met_child)
    print("Methylation matrix saved")
    
    location=[]
    for probe in probes:
        loc=map_chr[np.where(map_chr[:,3]==probe),1][0][0]
        location.append(loc)
        
    with open(out_dir+"probes_chr_"+chrom+".csv", "w") as f:
         writer = csv.writer(f)
         writer.writerows(zip(probes,location))
    
    print("Probes info saved")
