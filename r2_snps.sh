#!/bin/bash
#SBATCH --time=240:00:00
#SBATCH --mem=50Gb
#SBATCH -c 3
#SBATCH --array=1-21
#SBATCH --mail-type=ALL
#SBATCH --job-name=ld_snps
#SBATCH --output=logs/ld%a.log

unset SLURM_EXPORT_ENV

CHR=${SLURM_ARRAY_TASK_ID}

FULL_LIST=LD_fineMap/5matchSNPs/mQTLs_fineMap.txt
VCF_IN=LD_fineMap/2filterMendelian/vcf_filt/chr${CHR}.vcf
CHR_LIST=mQTLs_chr${CHR}.txt

awk -v chr=${CHR} '$1 == chr {print $2}' ${FULL_LIST} > ${CHR_LIST}

if [ ! -s ${CHR_LIST} ]; then
    echo "No SNPs found for chromosome ${CHR}, skipping"
    exit 0
fi

module load plink
# Remove duplicate RSIDs from VCF before running plink
plink2 --vcf ${VCF_IN} \
      --snps-only just-acgt \
      --rm-dup exclude-all \
      --make-bed \
      --out tmp_${CHR}

module load plink/1.90

# Then run LD on the deduplicated bed files
plink --bfile tmp_${CHR} \
      --ld-snp-list ${CHR_LIST} \
      --ld-window-r2 0.5 \
      --ld-window 999999 \
      --ld-window-kb 1000 \
      --r2 \
      --out ld_files/ld_${CHR}

# Clean up temp files
rm tmp_${CHR}.bed tmp_${CHR}.bim tmp_${CHR}.fam tmp_${CHR}.log
plink --vcf ${VCF_IN} --ld-snp-list ${CHR_LIST} --ld-window-r2 0.5 --ld-window 999999 --r2 --out ld_files/ld_${CHR}


