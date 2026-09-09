#!/bin/bash
module load bedtools2/2.30.0
#make fasta compatible; homo_sapiens_ancestor_GRCh38 comes from 
#https://ftp.ensembl.org/pub/release-106/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh38.tar.gz
for chr in {1..22}; do
  input="homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_${chr}.fa"
  output="homo_sapiens_ancestor_GRCh38/ancestral_${chr}.fa"

  sed -E "s/>ANCESTOR_for_chromosome:GRCh38:([0-9]+):([0-9]+):([0-9]+):1/>chr\1/" "$input" > "$output"
done

#use getfasta to extract ancestral seq
mkdir -p ancestral_SNPs_outputs_1KG_base_comp
for chr in {1..22}; do
  # Define paths
  fasta="homo_sapiens_ancestor_GRCh38/ancestral_${chr}.fa"
  coords="1KG_incl_singl_coordinates_base_comp/chr${chr}.bed"
  output="ancestral_SNPs_outputs_1KG_base_comp/ancestral_${chr}.fa.out"

  # Use bedtools to extract ancestral alleles
  bedtools getfasta -fi "$fasta" -bed "$coords" -fo "$output"
done

#process fasta output into .bed
mkdir -p bed_files_1KG_base_comp
for chr in {1..22}; do
  # Input .fa.out file
  fa_out="ancestral_SNPs_outputs_1KG_base_comp/ancestral_${chr}.fa.out"
  bed_output="bed_files_1KG_base_comp/ancestral_${chr}.bed"
  # Convert to BED format
  awk 'BEGIN {FS="[:-]"; OFS="\t"} 
     /^>/ {chrom=$1; start=$2; end=$3; next} 
     {print substr(chrom, 2), start-1, end, toupper($1)}' "$fa_out" > "$bed_output"
done

