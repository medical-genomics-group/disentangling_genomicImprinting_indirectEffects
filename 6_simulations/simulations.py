import numpy as np
import zarr
import os
import sys
import argparse


def index_r2(r2,limit1=0.2,limit2=1): #give back index of SNPs with certain r2 
    filt=r2[np.where((r2[:,2]>limit1) & (r2[:,2]<limit2))]
    r2_filt=filt[:,0].tolist() #first snp
    r2_filt+=filt[:,1].tolist() #second snp
    r2_filt=np.unique(r2_filt)
    return([int(x) for x in r2_filt])

def getSNP(noCausal,noSNP,freq=None,r2=None): #get list of indexes for SNPs with certain minor allele frequency and r2 to inject beta
    if freq and r2:
        freq_r2=list(set(freq).intersection(r2))
        snp=np.random.choice(freq_r2,noCausal,False)
        snp=np.take(noSNP,snp)
    elif freq:
        snp=np.random.choice(freq,noCausal,False)
        snp=np.take(noSNP,snp)
    elif r2:
        snp=np.random.choice(r2,noCausal,False)
        snp=np.take(noSNP,snp)
    else: #random
        snp=np.random.choice(noSNP,noCausal,False)
    return(snp.tolist())

def get_covMatrix(pType): #choose pattern of phenotype (maternal,paternal,complex)
    pType=np.abs(pType)
    cov=[]
    for i in pType:
        cov_i=[]
        for j in pType:
            cov_i.append(i*j)
        cov.append(cov_i)
    return(np.array(cov))

def getBetas(size,noCausal,h2,snps,pType=[1,0,0,0]): 
    beta=np.zeros(size)
    cov=get_covMatrix(pType)
    V=(h2/noCausal)*cov
    values=np.random.multivariate_normal([0,0,0,0],V,noCausal) #get coefficients from distribution
    values=values*pType #keep sign of signals (important for CI)
    var=np.var(values,axis=0) #actual variance of injected signals
    
    for snp,row in zip(snps,values):
        beta[snp]=row[0] #C
        beta[snp+1]=row[1]#M
        beta[snp+2]=row[2]#F
        beta[snp+3]=row[3]#I
  
    return (beta,var)


def main(pathX,path_freq,path_r2,path_snps,pType,h2,freq_no,noCausal,r2_no,fix,out_dir,out_name):
    group=4 #number of parameters
    ###### r2 load files and change snps names for index in snps_names
    r2=np.loadtxt(path_r2,dtype=object)
    r2=np.delete(r2,[0,1,3,4],axis=1)

    with open(path_snps, 'r') as file:
        snps_names= [line.strip() for line in file]

    delete=[] #delete snps that aren't in list
    for i in range(0,len(r2)):
        snp1=r2[i,0]
        snp2=r2[i,1]
        if not(snp1 in snps_names) or not(snp2 in snps_names):
            delete.append(i)

    r2=np.delete(r2,delete,axis=0)

    r2[:,0]=[snps_names.index(x) for x in r2[:,0]] #change names to index
    r2[:,1]=[snps_names.index(x) for x in r2[:,1]]

    ####
    match r2_no:
        case 0:#<0.2
            r2_ind=[(x in r2[:,0]) or (x in r2[:,1]) for x in range(0,len(snps_names))]
            r2_ind= [i for i, val in enumerate(r2_ind) if not val] #not related, min value 0.2
        case 1:#0.2-0.5
            r2_ind=index_r2(np.array(r2,dtype=float),0.2,0.5)
        case 2:#0.5-0.75
            r2_ind=index_r2(np.array(r2,dtype=float),0.5,0.75)
        case 3:#0.75-1
            r2_ind=index_r2(np.array(r2,dtype=float),0.75)
        case 4:#random
            r2_ind=None
        case _:
            r2_ind=None
            print("Warning: r2_no not identified random selection")
            
    ####################
    X=zarr.load(pathX)
    #STANDARDIZE X
    Xstd=(X-np.mean(X,axis=0))/(np.std(X,axis=0)) 

    noSNP=list(range(0,X.shape[1],group)) #starting index of snps C component

    #load MAF file
    freq=np.loadtxt(path_freq,skiprows=1,delimiter="\t",dtype=object) #in order with snps
    freq=np.delete(freq,[0,1,2,3,5],axis=1)
    freq=freq[:,0]
    freq=[float(x.split(":")[1]) for x in freq]
    limits=[0.05,0.1,0.3,0.5]

    ##########################################################################
    ### select SNPs and inject beta
    if freq_no==3: #random
       SNPs=getSNP(noCausal,noSNP,freq=None,r2=r2_ind)
       beta,V=getBetas(X.shape[1],noCausal,h2,SNPs,pType)
    else:
        min1=limits[freq_no]
        min2=limits[freq_no+1]
        max1=1-min1
        max2=1-min2
        
        ind_freq=[j for j in range(len(freq)) if (freq[j]>=min1 and freq[j]<=min2) 
             or (freq[j]<=max1 and freq[j]>=max2)] #indexes of snps within the bounces of the MAF
        
        SNPs=getSNP(noCausal,noSNP,ind_freq,r2_ind)
        MAF_selected=[freq[int(x/4)] for x in SNPs]
        beta,V=getBetas(X.shape[1],noCausal,h2,SNPs,pType)
        

    #inject signals in beta, calculate y
    g=np.matmul(Xstd,beta.T)

    #fix h2=V
    if fix:
        sigma2 = np.var(g) * ((1-h2)/h2)
    else:
        sigma2=np.sqrt(1-h2)
        
    e = np.random.normal(0, np.sqrt(sigma2), X.shape[0])
    y=g+e   

    print("Calculated h2 after adding noise:")
    h2_calculated=np.var(Xstd@beta)/np.var(y)
    print(h2_calculated)

    print("Calculated variance of injected signals:")
    print(V)

    #standardize y for models
    ystd=(y-np.mean(y))/np.std(y)

    file_name="_".join( str(x) for x in out_name )
    zarr.save(out_dir+"beta_"+file_name+'.zarr',beta)
    zarr.save(out_dir+"y_"+file_name+'.zarr',ystd)

    #SAVE Y AS TXT FOR JODIE
    np.savetxt(out_dir+"y_"+file_name+'.txt',ystd)

    #save info on snps selected
    #MAF AND NAME
    np.savetxt(out_dir+"Variance_"+file_name+'.txt',V) 
    np.savetxt(out_dir+"MAF_"+file_name+'.txt',MAF_selected) 
    
    print("beta done!")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Simulate CMFI signals on top of vcf chromosome file, outputs coefficients')
    parser.add_argument('--pathX', type=str, help='Path to zarr file of X matrix for chormosome',required=True)
    parser.add_argument('--path_freq', type=str, help='File .frq with MAF from SNPs',required=True)
    parser.add_argument('--path_r2', type=str, help='.vcor file with r2',required=True)
    parser.add_argument('--path_snps', type=str, help='File .txt with SNPs names to match r2 and MAF',required=True)
    
    parser.add_argument('--pType',type=list,default=[1,0,0,0],help='Simulated pattern [C,M,F,I]')
    parser.add_argument('--h2', type=float,default=0.2,help='h2 or V')
    parser.add_argument('--freq_no', type=int,default=2,help='MAF of selected SNPs 0=5%-10% 1=10%-30% 2=30%-50% 3=random default=2')
    parser.add_argument('--noCausal', type=int,default=5,help='Number of causal SNPs with signal !=0')
    parser.add_argument('--r2_no', type=int,default=0,help='r2 between selected SNPs 0=r2<0.2 1=[0.2,0.5] 2=[0.5,0.75] 3=[0.75,1] 4=random default=0')
    parser.add_argument('--fix', type=bool,default=True,help='Fix h2=V')

    
    parser.add_argument('--out_dir', type=str,required = True)
    parser.add_argument('--out_name', type=str,required = True) #out_name,pType,freq_no,noCausal,r2_no

    args = parser.parse_args()
    print(args)

    main(args.pathX,args.path_freq,args.path_r2,args.path_snps,args.pType,args.h2,args.freq_no,
         args.noCausal,args.r2_no,args.fix,args.out_dir,args.out_name)
    print("Done")
