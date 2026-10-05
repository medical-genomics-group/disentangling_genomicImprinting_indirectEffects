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

"""
Get candidates only indirect effects
"""
normVal=3.474051424481396
var_thr=0.05

path="screen_full"
screen2_dir="POE_mQTLs/outputs/screen2CMFI.csv"

path_snps="screen_full/X_matrix/out_x/snps_"

screen2=pd.read_csv(screen2_dir) 

filt=screen2[(screen2[["var varC","var varM","var varF","var varI"]]<normVal).all(axis=1)] #filtered by var of vars

#filter screen2 to those with significant varF or varM
m_probes=filt[filt["varM"]>var_thr]
f_probes=filt[filt["varF"]>var_thr]

## add patterns
m_probes.to_csv("POE_mQTLs/outputs/mat_probes.csv",index=False)
f_probes.to_csv("POE_mQTLs/outputs/pat_probes.csv",index=False)
    








