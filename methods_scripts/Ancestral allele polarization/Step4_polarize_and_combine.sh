#!/bin/bash
# polarize_and_combine.sh: using ancestral alleles, polarize mutations and
#re-calculate DAF

INDIR="bed_files_with_ancestral"
OUT="/path/to/output/all_chrs_polarized.bed"

# Header
echo -e "CHROM\tSTART\tEND\tMUTATION\tDAF_AFR\tDAF_EUR\tDAF_AMR\tDAF_EAS\tDAF_SAS\tDAF_GLOBAL" > ${OUT}

for CHR in {1..22}; do
    echo "Processing chr${CHR}..."
    zcat ${INDIR}/chr${CHR}_snps_ancestral.bed.gz \
    | awk 'BEGIN{OFS="\t"} {
        chrom      = $1
        start      = $2
        end        = $3
        ref        = $4
        alt        = $5
        daf_afr    = $6
        daf_eur    = $7
        daf_amr    = $8
        daf_eas    = $9
        daf_sas    = $10
        daf_global = $11
        anc        = $12

        # Skip missing/ambiguous ancestral
        if (anc == "." || anc == "N" || anc == "") next

        if (anc == ref) {
            # ALT is derived — DAF already reflects derived allele frequency
            mutation = ref">"alt
            d_afr    = daf_afr
            d_eur    = daf_eur
            d_amr    = daf_amr
            d_eas    = daf_eas
            d_sas    = daf_sas
            d_global = daf_global

        } else if (anc == alt) {
            # REF is derived — DAF reflects ancestral freq, so flip it
            mutation = alt">"ref
            d_afr    = (daf_afr    == "NA") ? "NA" : 1 - daf_afr
            d_eur    = (daf_eur    == "NA") ? "NA" : 1 - daf_eur
            d_amr    = (daf_amr    == "NA") ? "NA" : 1 - daf_amr
            d_eas    = (daf_eas    == "NA") ? "NA" : 1 - daf_eas
            d_sas    = (daf_sas    == "NA") ? "NA" : 1 - daf_sas
            d_global = (daf_global == "NA") ? "NA" : 1 - daf_global

        } else {
            next    # ancestral matches neither allele — skip
        }

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
            chrom, start, end, mutation,
            d_afr, d_eur, d_amr, d_eas, d_sas, d_global

    }' >> ${OUT}
done

echo "Done: $(wc -l < ${OUT}) lines (including header)"
echo "Output: ${OUT}"