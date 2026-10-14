#!/bin/bash
#SBATCH --time=240:00:00
#SBATCH --mem=50Gb
#SBATCH -c 3
#SBATCH --array=1-400
#SBATCH --mail-type=ALL
#SBATCH --job-name=sim
#SBATCH --output=logs/sim_%a.log

unset SLURM_EXPORT_ENV
id=${SLURM_ARRAY_TASK_ID}

module load R
module load python/3.11.4
module load openmpi/4.1.4

#patterns
patterns=([-1,0,0,1] [1,0,0,-1] [-1,0,0,-1] [1,0,0,1] [0,0,0,1] [0,0,0,-1] [1,1,0,0] [1,0,1,0])
names_pat=("M-" "M+" "F-" "F+" "C+" "C-" "CM" "CF")

index=$(( ($id - 1) % ${#patterns[@]} ))
PAT=${patterns[$index]}
PAT_name=${names_pat[$index]}

echo "ID "$id
echo "pattern: "$PAT

mkdir -p out

mkdir -p out_beta_y
mkdir -p out_lasso
mkdir -p out_jodie
mkdir -p out_susieR

module load python/3.12

#########################################
beta_dir="simulation.py"
x_dir="sim_whole22_parametersLasso/chr22_X/x_22.zarr/"
freq_dir="sim_whole22_parametersLasso/chr22_X/freq_22.frq"
r2_dir="sim_whole22_parametersLasso/chr22_X/r2_22.vcor"
snps_dir="sim_whole22_parametersLasso/chr22_X/snps22.txt"


jodie="jodie/JODIE-main/"
susie="susieR/"
lasso="lasso/run_lasso_parameters.py"
path=$(pwd)

######################
#1.beta
name=${id}_${PAT_name}
srun python ${beta_dir} --pathX ${x_dir} --path_freq ${freq_dir} --path_r2 ${r2_dir}\
 --path_snps ${snps_dir} --ptype ${PAT} --out_dir out_beta_y/ --out_name ${name} 
echo 'Betas finished'


####################
#2.run lasso
srun python ${lasso} ${x_dir} out_beta_y/y_${name}.zarr out_lasso/ ${name}
echo 'lasso finished'

####################
#3.run jodie
# 3.0.calc_xtx.py #just needs to run once
#srun python ${jodie}calc_xtx.py --n 2438 --p 16466 --k 4 --pheno pheno --x ${x_dir} --dir out_jodie/ 
#echo 'calc_xtx jodie finished'

#alredy calculated
xtx="X_matrix/outXTX/XtX_22.zarr/"

# 3.1.parental_gibbs_sampler.py
mkdir -p out_jodie/${name}
#mod_parental
srun python ${jodie}jodie.py --iters 1000 --burnin 500 --n 2438 --p 16466 --k 4 --g 16466 --x ${x_dir} --xtx ${xtx} --y out_beta_y/y_${name}.txt --resultdir out_jodie/${name}/

echo 'gibbs sampler jodie finished'

echo 'jodie finished'

###############################

#4. susie
srun Rscript ${susie}susieR.R ${x_dir} ${path}/out_beta_y/ ${name} ${path}/out_susieR/

echo "susie finished"
