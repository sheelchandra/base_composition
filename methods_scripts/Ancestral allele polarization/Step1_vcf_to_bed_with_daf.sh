#!/bin/bash
# vcf_to_bed_with_daf.sh
# For each chr, create a BED file with REF, ALT, and per-superpopulation DAF that can be used
#later for polarization and for calculating GC%
# Output: bed_files/chr{N}_snps.bed.gz

#BSUB -J vcf2bed[1-22]
#BSUB -o logs/vcf2bed_%I.out
#BSUB -e logs/vcf2bed_%I.err
#BSUB -n 3
#BSUB -M 160000
#BSUB -R "rusage[mem=160000]"
module load python
module load bcftools
module load htslib
module load samtools

set -euo pipefail

CHR=${LSB_JOBINDEX}
VCF="all/chr${CHR}_all.vcf.gz" #putatively neutral VCFs for each chr (no MAF filter)
SAMPLE_INFO="non_admixed_1KG_sample_IDs.txt"
OUTDIR="bed_files"
mkdir -p ${OUTDIR}
OUT="${OUTDIR}/chr${CHR}_snps.bed"

echo "[$(date)] Processing chr${CHR}..."

python - <<EOF
import gzip, sys
import numpy as np
from collections import defaultdict

vcf      = "${VCF}"
sampfile = "${SAMPLE_INFO}"
out      = "${OUT}"

POPS = ["AFR", "EUR", "AMR", "EAS", "SAS"]

sample_pop = {}
with open(sampfile) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 3:
            sample_pop[parts[0]] = parts[2]

with gzip.open(vcf, "rt") as f, open(out, "w") as o:

    o.write("\t".join(["#CHROM", "START", "END", "REF", "ALT"] +
                      [f"DAF_{p}" for p in POPS] +
                      ["DAF_GLOBAL"]) + "\n")

    for line in f:
        if line.startswith("##"):
            continue

        if line.startswith("#CHROM"):
            all_samples = line.strip().split("\t")[9:]
            pop_indices = defaultdict(list)
            for i, s in enumerate(all_samples):
                pop = sample_pop.get(s)
                if pop in POPS:
                    pop_indices[pop].append(i)
            continue

        fields = line.strip().split("\t")
        chrom  = fields[0]
        pos    = int(fields[1])
        ref    = fields[3]
        alt    = fields[4]

        if len(ref) != 1 or len(alt) != 1:
            continue
        if ref not in "ACGT" or alt not in "ACGT":
            continue

        genotypes = fields[9:]

        # Per-population DAF
        pop_dafs = []
        for pop in POPS:
            derived = 0
            total   = 0
            for i in pop_indices[pop]:
                gt = genotypes[i].split(":")[0].replace("|", "/")
                if "." in gt:
                    continue
                for a in gt.split("/"):
                    if a == "1":
                        derived += 1
                    total += 1
            daf = derived / total if total > 0 else "NA"
            pop_dafs.append(f"{daf:.6f}" if daf != "NA" else "NA")

        # Global DAF — all non-admixed individuals combined
        global_derived = 0
        global_total   = 0
        for pop in POPS:
            for i in pop_indices[pop]:
                gt = genotypes[i].split(":")[0].replace("|", "/")
                if "." in gt:
                    continue
                for a in gt.split("/"):
                    if a == "1":
                        global_derived += 1
                    global_total += 1

        global_daf = f"{global_derived / global_total:.6f}" if global_total > 0 else "NA"

        o.write("\t".join([chrom, str(pos - 1), str(pos), ref, alt] +
                          pop_dafs +
                          [global_daf]) + "\n")

print(f"Done: {out}", file=sys.stderr)
EOF

# Compress and index
bgzip ${OUT}
tabix -p bed ${OUT}.gz

echo "[$(date)] Done chr${CHR} -> ${OUTDIR}/chr${CHR}_snps.bed.gz"
