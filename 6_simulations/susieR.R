`%ni%` = Negate(`%in%`)
library(argparse)
N_CORE=32

library(susieR)
library(Rarr)
##########################

#1. Load X and Y files
args=commandArgs(trailingOnly=T)
xdir=args[1]
path=args[2]
name=args[3]
out_dir=args[4]


y_dir=paste(path,"y_",name,".zarr",sep="")

x<-read_zarr_array(xdir)
y<-read_zarr_array(y_dir)

beta_name=paste("beta_",name,".zarr",sep="")
pip_name=paste("pip_",name,".zarr",sep="")

dir_beta_zarr=file.path(out_dir,beta_name)
dir_pip_zarr=file.path(out_dir,pip_name)

# Solve 
res <- susie(x,y,L=10) #Ls maximum number of non-zero effects
betas= coef(res)
betas=array(betas,dim=c(length(betas),1))

print("Coefficients calculated")

#save calculated betas
print(dir_beta_zarr)
write_zarr_array(x = betas, zarr_array_path=dir_beta_zarr, chunk_dim=dim(betas))


pip=res$pip
pip=array(pip,dim=c(length(pip),1))
write_zarr_array(x=pip, zarr_array_path=dir_pip_zarr, chunk_dim=dim(pip))

print("Files saved")

