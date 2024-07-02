#!/bin/bash

# Get assembly statistics for SaFo, lake trout (SaNa), rainbow trout (OMyk), brown trout (SaTr), Atlantic salmon (SaSa), whitefish (CoCl) and Dolly Varden (SaMa) close species
# Run in a conda env where assembly-stats has been installed

# srun -c 1 -p small --mem=20G --time=1-00:00:00 -J 03_species_comp -o log/03_species_comp_%j.log /bin/sh ./01_scripts/03_species_comp.sh &

# VARIABLES

COMP_DIR="species_comparison/GCF"

SAFO="$COMP_DIR/GCF_029448725.1_ASM2944872v1_genomic.fna"
SANA="$COMP_DIR/GCF_016432855.1_SaNama_1.0_genomic.fna"
OMYK="$COMP_DIR/GCF_013265735.2_USDA_OmykA_1.1_genomic.fna"
SATR="$COMP_DIR/GCF_901001165.1_fSalTru1.1_genomic.fna"
SASA="$COMP_DIR/GCF_905237065.1_Ssal_v3.1_genomic.fna"
COCL="$COMP_DIR/GCF_020615455.1_ASM2061545v1_genomic.fna"
SAMA="$COMP_DIR/GCF_002910315.2_ASM291031v2_genomic.fna"

SASA_IC="$COMP_DIR/GCF_000233375.1_ICSASG_v2_genomic.fna"   #ICSASG assembly v2

# LOAD REQUIRED MODULES
module load samtools/1.15
module load bedtools/2.30.0
module load BBTools/36.92


# Run assembly-stats
assembly-stats -t $SAFO $SANA $OMYK $SATR $SASA $COCL $SAMA $SASA_IC


# Loop over assemblies
for i in $SAFO $SANA $OMYK $SATR $SASA $COCL $SAMA $SASA_IC; 
do
  # Index genome
  echo $i
  samtools faidx $i
  
  # Run stats.sh from BBTools
  stats.sh $i > "${i%.*}".BBTools.stats.txt
  
  # Get number of bp anchored into chromosomes
  ## Make bed and extract chr: Unplaced contigs names start with scaf
  less "$i".fai | awk -F '\t' '{printf("%s\t0\t%s\n",$1,$2);}' | grep -E "^CM|HG|LR|NC" > "${i%.*}".chrs.bed
  less "${i%.*}".chrs.bed | wc -l 
  
  ## Remove unplaced scaffolds from fasta
  bedtools getfasta -fi "$i" -bed "${i%.*}".chrs.bed -fullHeader > "${i%.*}".chrs.fasta
  
  ## Count bases in CHR
  CHR_BP=$(less "${i%.*}".chrs.fasta | grep -E '^[AaCcTtGgNn]+$' | perl -pe 's/[[:space:]]//g' | wc -c)
  echo "$CHR_BP bases in chromosomes"
  
  ## Count total bases in assembly
  GENOME_BP=$(less $i | grep -E '^[AaCcTtGgNn]+$' | perl -pe 's/[[:space:]]//g' | wc -c)
  echo "$GENOME_BP total base pairs in fasta"
  
  ## Divide both
  echo "$CHR_BP / $GENOME_BP"
done