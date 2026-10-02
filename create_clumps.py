import numpy as np
import pandas as pd
import os

path="LD_fineMap/5matchSNPs/"

ld_df=[]

for chro in range(0,22):
    chro=str(chro)
    ld_file=path+"ld_files/"+"ld_"+str(chro)+".ld"
    if os.path.exists(ld_file):
       ld=pd.read_csv(ld_file, sep=r"\s+",usecols=["SNP_A","SNP_B","R2"])
       ld_df.append(ld)
    else:
       print("ld does not exist for chromosome ",chro)

ld_df= pd.concat(ld_df,ignore_index=True)
ld_df.to_csv(path+"r2_mQTL_fineMap.csv",index=False)
       
