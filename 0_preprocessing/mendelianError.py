import pandas as pd
import numpy as np
from pandas_plink import read_plink1_bin
import zarr
import sys
import os
import math
from collections import Counter
from make_x.py import getGeno, getGeno_Haplo
import argparse

def main(chro,path_trios,path_ids,path_genotype,X,out_dir):
    # load trios and genotype
    vcf=np.genfromtxt(path_genotype+".vcf", delimiter='\t', comments='#', dtype=object)
    vcf= np.delete(vcf,[0, 1, 3, 4, 5, 6, 7, 8], axis=1)  # Keep only needed columns

    ids=np.loadtxt(path_ids,delimiter="\t",dtype=int)
    ids=ids.tolist()

    trios=np.loadtxt(path_trios,dtype=int)

    snps=vcf[:,0]
    vcf=np.delete(vcf,0,1) #vcf[snp,individual]

    #load ids
    print('Files loaded')

    error=[]
    error_phasing=[]
    for i in range(0,len(trios)): #per child
        for j in range(0,len(vcf)): #per snp
            loc_child=ids.index(trios[i,0])
            loc_mother=ids.index(trios[i,2])
            loc_father=ids.index(trios[i,1]) #father
        
            
            child=getGeno(vcf[j,loc_child].decode(),True) #child, mother, father
            mother=getGeno(vcf[j,loc_mother].decode())
            if X==True:
                father=getGeno(vcf[j,loc_father].decode())
            else:
                father=getGeno(vcf[j,loc_father].decode(),False)
            parental=[mother,father]
            parental.sort()
            
            snp=snps[j].decode()
            
            match parental:
                case [0]: #aa
                    if child[0]!=0:
                        error.append(snp)
                case [0,1]:
                    if child[0]==2: #cannot be AA
                        error.append(snp)
                case [0,2]:
                    if child[0]!=1:
                        error.append(snp)
                case [1,2]:
                    if child[0]==0: #cannot be aa
                        error.append(snp)
                case [2]:
                    if child[0]!=2:
                        error.append(snp)
                        
            #mother
            if mother==2 and child[1]==0:
                error_phasing.append(snp)
            if mother==0 and child[1]==1:
                error_phasing.append(snp)
                
            #father
            if father==2 and child[2]==0:
                error_phasing.append(snp)
            if father==0 and child[2]==1:
                error_phasing.append(snp)
                
                

    print("Finished with errors")

    total_error=len(error)
    print("Total number of errors:",total_error)
    percent=(total_error*100)/(len(trios)*len(vcf))
    print("Percentage of errors: ",percent)

    with open(out_dir+"error_all.txt", "a+") as f:
          s = ','.join(str(x) for x in [chr,total_error,percent])
          f.write(s+"\n")

    #####
    tests=len(trios)
    threshold=math.ceil(tests/100) #1% of error threshold 

    most_common_errors={}
    res=Counter(error)
    most_common_errors[chro]=res.most_common(1500)

    error_array=np.array(res.most_common(1500))
    snps_freq=[int(x) for x in error_array[:,1] if int(x)>=threshold]
    to_filter=list(range(0,len(snps_freq))) #index
    to_filter=error_array[to_filter,0]
        
    with open(out_dir+"snps_filter_"+chro+".txt", mode='wt', encoding='utf-8') as myfile:
         myfile.write('\n'.join(to_filter))
         
    ### ERROR PHASING
    total_error=len(error_phasing)
    print("Total number of errors:",total_error)
    percent=(total_error*100)/(len(trios)*len(vcf)*2)
    print("Percentage of errors: ",percent)


    tests=len(trios)
    threshold=math.ceil(tests/100) #1% of error threshold

    most_common_errors={}
    res=Counter(error_phasing)
    most_common_errors[chro]=res.most_common(1500)

    error_array=np.array(res.most_common(1500))
    snps_freq=[int(x) for x in error_array[:,1] if int(x)>=threshold]
    to_filter=list(range(0,len(snps_freq))) #index
    to_filter=error_array[to_filter,0]
        
    with open(out_dir+"snps_filter_phasing_"+chro+".txt", mode='wt', encoding='utf-8') as myfile:
         myfile.write('\n'.join(to_filter))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get mendelian error including phasing')
    parser.add_argument('--chro', type=str, default="1",help='Chromosome',)
    parser.add_argument('--path_trios', type=str, help='Dir to trios.ped', required = True)
    parser.add_argument('--path_ids', type=str, help='Dir to ids of vcf', required = True)
    parser.add_argument('--path_genotype', type=str, help="Dir to vcf of chromosome",required = True)
    parser.add_argument('--X', type=bool, default=False,help="Is it the X chromosome?/Haploid for fathers")
    parser.add_argument('--out_dir', type=str, default="", help="Dir for outputs")
    
    args = parser.parse_args()
    print(args)

    main(chro=args.chro,path_trios=args.path_trios,path_ids=args.path_ids,path_genotype=args.path_genotype,
         X=args.X,out_dir=args.out_dir)
    print("Done")
