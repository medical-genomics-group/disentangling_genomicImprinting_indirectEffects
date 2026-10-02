"""
Created on Wed Jun  4 11:45:26 2025

@author: madamega
make trios, duos and siblings lists
"""

import numpy as np
import pandas as pd
import sys

pd.options.mode.chained_assignment = None 

def locatePOID(PO,id1):
    withID=kin[np.where((kin[:,0]==id1) | (kin[:,1]==id1))[0]]
    return(withID)


#output path
out_dir="out_family/"

### load king, variables and linkID files:
kin_path="0_grouping/GS_king.kin"
kin=pd.read_csv(kin_path,delimiter="\t")
kin=kin.drop(["FID","N_SNP","Z0","Phi","HetHet","HetConc","HomIBS0","IBD2Seg","IBD1Seg","PropIBD","Error"],axis=1)

idlink_path='0_grouping/GS_idlink.csv'
idLink=pd.read_csv(idlink_path,dtype=object)

covars=pd.read_csv("covars.csv") #file with covariances 
idLink=pd.read_csv("1.matrix/Y/idlink.csv")
covars=covars.rename(columns={"Sample_Sentrix_ID": "X"})
covars=covars.merge(idLink,on='X')

### merge sex and age to kin df
kin=pd.merge(kin,covars[["Sample_Name","age","sex"]],left_on="ID1",right_on="Sample_Name",how="left")
kin=kin.rename(columns={"age":"ID1_Age","sex":"ID1_sex"})
kin=pd.merge(kin,covars[["Sample_Name","age","sex"]],left_on="ID2",right_on="Sample_Name",how="left")
kin=kin.rename(columns={"age":"ID2_Age","sex":"ID2_sex"})

############################################################### ONLY CONSIDERING SAMPLES WITH SEX AND AGE AVAILABLE!!!!!
#drop nans
kin=kin.dropna()

#save this samples
kin.to_csv(out_dir+"kin_covars.csv")
print("Saved kin + covars, samples considered data frame")

###################################################################################
### trios and duos
##set thresholds
po_down=1/(2**(5/2))
po_up=1/(2**(3/2))
cut_IBS0=0.0012
cut_age=15

## select only parent-offspling pairs
PO=kin[(kin['Kinship'] >= po_down) & (kin['Kinship']<po_up) & (kin['IBS0']<cut_IBS0)]
print("sanity check unique types of infered type from king based on PO filters:")
print(np.unique(PO["InfType"]))

### RULES FOR PARENTS
#kinship>po_down, kinship<po_up, IBS0<cut_IBS0
print("filter kinship>po_down: ",po_down)
print("filter kinship>po_up: ",po_up)
print("filter IBS0<cut_IBS0: ",cut_IBS0)

## diff in age between parent-offspring is more than 15 years
## make sure parents dont have the same genetic sex
trios=pd.DataFrame({"child":[],"M":[],"F":[]}) #CHILDREN, MALE, FEMALE !!!

PO["diff"]=PO["ID1_Age"]-PO["ID2_Age"]

dropped_PO=PO[(np.abs(PO["diff"])<cut_age)]
print("Number of entries dropped differences of ages <15: ",len(dropped_PO))
PO=PO[(np.abs(PO["diff"])>=cut_age)] # filter difference in age <15

sameSex=[]
for ind,row in PO.iterrows():
    if row["diff"]<0: #id2 is parent
        child=str(row["ID1"])
        parent=str(row["ID2"])
        sex_p=row["ID2_sex"]
    else: #id1 is parent
        child=str(row["ID2"])
        parent=str(row["ID1"])
        sex_p=row["ID1_sex"]
       
    #children already in df trios?
    ind_t=trios.index[trios["child"]==child].tolist()
    
    if ind_t: #child already in trios as child
        #check if parent of same sex already there
        if len(ind_t)>1:
            print("Error, child appears many times")
        else:
            ind_t=ind_t[0]
        if trios.iloc[ind_t][sex_p].isnumeric(): #parent with that sex already there
            print("Two parents with the same sex")
            print("index: ",ind_t,"child: ",child)
            sameSex.append(child)
            #if difference in age is exactly the same, they are the same parent #add in new version
        else:
            trios.at[ind_t,sex_p]=parent
        
    else: #child not in trios df
        if sex_p=="M":
            ind_p=1
        else:
            ind_p=2
        new_trio=[child,"None","None"]
        new_trio[ind_p]=parent
        
        trios.loc[len(trios)]=new_trio

########
#two parents same sex
for i in sameSex:   
    i=int(i)
    a=PO.index[PO["ID1"]==i].tolist()
    b=PO.index[PO["ID2"]==i].tolist()
    a=a+b
    filt_PO=PO.loc[a]
    print(filt_PO)

###############################################################33
###OUTPUT
#save dataframe
trios.to_csv(out_dir+"trios_all.csv")

### trios.ped
## childID  fatherID motherID

trios.replace("None", np.nan, inplace=True)
complete_trios=trios.dropna()

#save trios
complete_trios.to_csv(out_dir+"trios.ped",sep="\t",header=False,index=False)
complete_trios=complete_trios.astype(int)

#ids
list_trios=complete_trios["child"].tolist()+complete_trios["M"].tolist()+complete_trios["F"].tolist()
list_trios=np.unique(list_trios)

with open(out_dir+'trios_ids.txt','w') as f:
    for line in list_trios:
        f.write(f"{line}\n")
        

### duos.ped
## childID fatherID/NAN motherID/NAN

duos=trios.loc[~trios.index.isin(trios.dropna().index)]
duos.replace(np.nan,"NA",inplace=True)
duos.to_csv(out_dir+"duos.ped",sep="\t",header=False,index=False)
duos.replace("NA",0,inplace=True)

duos=duos.astype(int)
list_duos=duos["child"].tolist()+duos["M"].tolist()+duos["F"].tolist()
list_duos=list(np.unique(list_duos))
list_duos.remove(0)

### siblings
sib=kin[(kin['Kinship'] >= po_down) & (kin['Kinship']<po_up) & (kin['IBS0']>cut_IBS0)]
print("sanity check unique types of infered type from king based on sib filters:")
print(np.unique(sib["InfType"]))

#siblings most have a difference in age lower than 15 years
sib["diff"]=np.abs(sib["ID1_Age"]-sib["ID2_Age"])

dropped_sib=sib[(sib["diff"])>15]
print("siblings dropped by difference in age >15")
print(len(dropped_sib))

sib=sib[(np.abs(sib["diff"])<=15)] # filter difference in age <15
print("sanity check unique types of infered type from king based on sib filters:")
print(np.unique(sib["InfType"]))
## ped file siblings
# ID1, ID2
#save df
sib.to_csv(out_dir+"siblings.csv")
ped_sib=sib[["ID1","ID2"]]
ped_sib.to_csv(out_dir+"siblings.ped",sep="\t",header=False,index=False)

#save unique list siblings
list_sib=ped_sib["ID1"].tolist()+ped_sib["ID2"].tolist()
list_sib=np.unique(list_sib)

with open(out_dir+"sib_list.txt", 'w') as fp:
    for item in list_sib:
        # write each item on a new line
        fp.write("%s\n" % item)

######
#extract siblings and duos not in trios
sib_not_trios=list(set(list_sib)-set(list_trios))
print("number samples in siblings and not in trios: ",len(sib_not_trios))
with open(out_dir+"sib_not_tios.txt", 'w') as fp:
    for item in sib_not_trios:
        # write each item on a new line
        fp.write("%s\n" % item)
        
duo_not_trios=list(set(list_duos)-set(list_trios))
with open(out_dir+"duo_not_trios.txt", 'w') as fp:
    for item in duo_not_trios:
        # write each item on a new line
        fp.write("%s\n" % item)


## dataframe sibs and duos not in trios
duos_notTrios=duos[~duos[['child','M','F']].isin(list_trios).any(axis=1)] #not in list_trios
sibs_notTrios=sib[~sib[['ID1','ID2']].isin(list_trios).any(axis=1)]

duos_notTrios.replace(0,"NA",inplace=True)

duos_notTrios.to_csv(out_dir+"duos_notTrios.ped",sep="\t",header=False,index=False)
sibs_notTrios[["ID1","ID2"]].to_csv(out_dir+"sibs_notTrios.ped",sep="\t",header=False,index=False)

### samples in covars not in trios
ids_covars=covars["Sample_Name"].tolist()
ids_notTrios=list(set(ids_covars)-set(list_trios))
ids_notTrios=pd.DataFrame(ids_notTrios,columns=["id"])

path_all="linearRegression_test/9proof/h2_siblings/SNPs/replicate/GS_GWAS_allIDS.txt"
all_ids=pd.read_csv(path_all,header=None,sep=" ",names=["famId","id"])

ids_notTrios=pd.merge(ids_notTrios,all_ids,on="id")
ids_notTrios=ids_notTrios.iloc[:,[1,0]]

ids_notTrios.to_csv(out_dir+"ids_not_trios.txt",sep="\t",index=False,header=None)


