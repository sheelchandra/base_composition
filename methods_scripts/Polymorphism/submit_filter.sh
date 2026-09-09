#!/bin/bash
# submit_filter.sh — submit as: bsub < submit_filter.sh

#BSUB -J filter_vcfs[1-2]       # array job, one per chromosome
#BSUB -o logs/filter_%I.out      # %I = array index
#BSUB -e logs/filter_%I.err
#BSUB -n 2                      
#BSUB -M 160000                    
#BSUB -R "rusage[mem=160000]"
module load bcftools
set -euo pipefail

mkdir -p logs all common rare
VCF_DIR="/path/to/VCFs/1KG_hg38/unphased_GT_vcfs"
SAMPLES="non_admixed_1KG_IDs.txt" #to only include unrelated, non-admixed individuals
REGIONS="regions_sorted.bed.gz" #pre-made file using RefSeq, 1KG Strictmask, and PhastCons coordinates
CHR=${LSB_JOBINDEX}              # LSF sets this to 1..22 automatically

INPUT=$VCF_DIR/"chr${CHR}.vcf.gz"
echo "[$(date)] Processing chr${CHR}..."

# Step 1: Filter → all/
bcftools view \
    --samples-file ${SAMPLES} \
    --force-samples \
    ${INPUT} \
| bcftools view \
    --min-alleles 2 --max-alleles 2 \
    --type snps \
| bcftools filter \
    --include 'F_MISSING < 0.2' \
| bcftools view \
    --targets-file ${REGIONS} \
| bcftools +fill-tags -- -t AF,MAF \
| bcftools view -Oz -o all/chr${CHR}_all.vcf.gz

bcftools index --tbi all/chr${CHR}_all.vcf.gz

# Step 2: Split common / rare (run in background in parallel within job)
bcftools view -i 'INFO/MAF >= 0.05' \
    all/chr${CHR}_all.vcf.gz \
    -Oz -o common/chr${CHR}_common.vcf.gz &

bcftools view -i 'INFO/MAF < 0.05' \
    all/chr${CHR}_all.vcf.gz \
    -Oz -o rare/chr${CHR}_rare.vcf.gz &

wait

bcftools index --tbi common/chr${CHR}_common.vcf.gz
bcftools index --tbi rare/chr${CHR}_rare.vcf.gz

echo "[$(date)] Done chr${CHR}"

