#!/bin/bash

# Get stats close species

# srun -c 1 -p small --mem=20G --time=1-00:00:00 -J species_comp -o log/species_comp_%j.log /bin/sh ./01_scripts/species_comp.sh &

# VARIABLES
FINAL_NCBI="08_final_NCBI/GCA_029448725.1_ASM2944872v1_genomic.fna"

COMP_DIR="species_comparison"
SANA="$COMP_DIR/GCA_016432855.1_SaNama_1.0_genomic.fna"
OMYK="$COMP_DIR/GCA_013265735.3_USDA_OmykA_1.1_genomic.fna"
SATR="$COMP_DIR/GCA_901001165.1_fSalTru1.1_genomic.fna"
SASA="$COMP_DIR/GCA_905237065.2_Ssal_v3.1_genomic.fna"
COCL="$COMP_DIR/GCA_020615455.1_ASM2061545v1_genomic.fna"


# LOAD REQUIRED MODULES
module load samtools/1.15
module load bedtools/2.30.0
module load BBTools

# First index genomes
#samtools faidx $SANA
#samtools faidx $OMYK
#samtools faidx $SATR
#samtools faidx $SASA
#samtools faidx $COCL

# Run assembly-stats
#assembly-stats -t $FINAL_NCBI $SANA $OMYK $SATR $SASA $COCL

# Run stats.sh
#stats.sh $SANA > "${SANA%.*}".stats.txt
#stats.sh $OMYK > "${OMYK%.*}".stats.txt
#stats.sh $SATR > "${SATR%.*}".stats.txt
#stats.sh $SASA > "${SASA%.*}".stats.txt
#stats.sh $COCL > "${COCL%.*}".stats.txt

# Get prop of bp in chr
## SANA
### Make bed and extract chr: Unplaced contigs names start with scaf
#less "$SANA".fai | awk -F '\t' '{printf("%s\t0\t%s\n",$1,$2);}' | grep -v 'JAE' > "${SANA%.*}".chrs.bed
### Remove unplaced scaffolds from fasta
#bedtools getfasta -fi "$SANA" -bed "${SANA%.*}".chrs.bed -fullHeader > "${SANA%.*}".chrs.fasta
### Count bases in CHR
#CHR_BP=$(less "${SANA%.*}".chrs.fasta | grep -E '^[AaCcTtGgNn]+$' | perl -pe 's/[[:space:]]//g' | wc -c)
#echo "$SANA"
#echo "$CHR_BP bases in chromosomes"
### Count bases
#GENOME_BP=$(less $FINAL_NCBI | grep -E '^[AaCcTtGgNn]+$' | perl -pe 's/[[:space:]]//g' | wc -c)
#echo "$GENOME_BP total base pairs in fasta"
#echo "$CHR_BP / $GENOM_BP"



for i in $SANA $OMYK $SATR $SASA $COCL; do
  echo $i
  ## Index
  samtools faidx $i
  ## Make bed and extract chr: Unplaced contigs names start with scaf
  less "$i".fai | awk -F '\t' '{printf("%s\t0\t%s\n",$1,$2);}' | grep -E "^CM|HG|LR" > "${i%.*}".chrs.bed
  less "${i%.*}".chrs.bed | wc -l 
  ## Remove unplaced scaffolds from fasta
  bedtools getfasta -fi "$i" -bed "${i%.*}".chrs.bed -fullHeader > "${i%.*}".chrs.fasta
  # Count bases in CHR
  CHR_BP=$(less "${i%.*}".chrs.fasta | grep -E '^[AaCcTtGgNn]+$' | perl -pe 's/[[:space:]]//g' | wc -c)
  echo "$CHR_BP bases in chromosomes"
  ## Count bases
  GENOME_BP=$(less $i | grep -E '^[AaCcTtGgNn]+$' | perl -pe 's/[[:space:]]//g' | wc -c)
  echo "$GENOME_BP total base pairs in fasta"
  echo "$CHR_BP / $GENOME_BP"
done