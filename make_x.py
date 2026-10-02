import numpy as np
import zarr
import sys
import os

chr = str(sys.argv[1])
out_dir = str(sys.argv[2])
pathInfo = str(sys.argv[3])
ids_file = str(sys.argv[4])
trios_file = str(sys.argv[5])
haploid= bool(sys.argv[6])

###############################################################################
#### Load VCF
info = np.genfromtxt(pathInfo + chr + ".vcf", delimiter='\t', comments='#', dtype=object)
#snps
snps=info[:,1:3]
info = np.delete(info, [0, 1, 3, 4, 5, 6, 7, 8], axis=1)  # Keep only needed columns

info=np.delete(info,0,1)
print('VCF file loaded')

###
deco_snps = np.vectorize(lambda x: x.decode('utf-8') if isinstance(x, bytes) else x)(snps)
np.savetxt(out_dir+"snps_"+chr+".txt", deco_snps, fmt='%s', delimiter='\t')

print("SNPs txt file saved: ",out_dir)

# Load IDs efficiently
ids = np.genfromtxt(ids_file, delimiter=" ", dtype=int)
ids_dict = {id_: idx for idx, id_ in enumerate(ids)}

# Load trio information
trios = np.loadtxt(trios_file, dtype=int)

print('Files loaded')

###############################################################################
#### Genotype
def getGeno_Haplo(snp_idx,sample):
    s_ind=ids_dict.get(sample,None)
    variant=info[snp_idx,s_ind].decode()
    x=variant.split(":")[0]
    return(x)

def getGeno(snp_idx, sample, child=False):
    s_ind = ids_dict.get(sample, None)  # Use dictionary lookup 
    variant = info[snp_idx, s_ind].decode()  # Convert bytes to string
    alleles = variant.split('|')
    # Extract allele information
    alleleF = int(alleles[0])
    alleleM = int(alleles[1].split(':')[0])
    x = alleleF + alleleM

    if child:
        if x==1:
            if alleleF==1:
                i=-1
            elif alleleM==1:
                i=1
        else:
            i=0
        return (x,i)
    else:
        return x


def xMatrix(matrix,no_snps,trio): 
    X=np.zeros([len(trio),no_snps*len(matrix)]) #shape=(n,j=p*k)
    for ind_c in range(0,len(trio)): #per child len(trio)
        ind=0
        for j in range(0,no_snps): #per snp nsnp
            m=getGeno(j,trio[ind_c,2])#mother 
            if haploid==True:
                f=getGeno_Haplo(j,trio[ind_c,1])
            else:
                f=getGeno(j,trio[ind_c,1])#father
            c,i=getGeno(j,trio[ind_c,0],True)
            geno=[c,m,f,i]                  

            for k in matrix: #CMFI
                X[ind_c,ind]=geno[k]
                ind+=1
    return X

###############################################################################
#### Execution
matrix=[0,1,2,3] #c,m,f,i,
print('Constructing X')
   
no_snps=len(snps)
X=xMatrix(matrix,no_snps,trios) #noComponents,n,n*s,ped file trios
print('X matrix done')

zarr.save(out_dir+"x_"+chr+'.zarr',X)
print("X matrix saved")

print('X matrix n: ',X.shape[0])
print('X matrix j variables: ',X.shape[1])
print("X matrix k groups: ",len(matrix),matrix)
print("X matrix p markers (snps):",no_snps)

