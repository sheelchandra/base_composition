#!/bin/bash
#BSUB -J coords[1-22]
#create coordinates directory for each SNP that can be used to extract ancestral alleles
set -euo pipefail

CHR=${LSB_JOBINDEX}
OUTDIR="1KG_incl_singl_coordinates_base_comp"
mkdir -p ${OUTDIR}
#bed files come from vcf_to_bed_daf.sh
zcat bed_files/chr${CHR}_snps.bed.gz \
    | grep -v "^#" \
    | awk '{print $1"\t"$2"\t"$3}' \
    > ${OUTDIR}/chr${CHR}.bed

echo "Done chr${CHR}: $(wc -l < ${OUTDIR}/chr${CHR}.bed) sites"

