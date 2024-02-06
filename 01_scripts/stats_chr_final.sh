#!/bin/bash

# Get stats for chromosomes

# srun -c 1 -p small --mem=20G --time=1-00:00:00 -J chr_stats_scaf2 -o log/chr_stats_scaf2_%j.log /bin/sh ./01_scripts/stats_chr_scaf2.sh &

# VARIABLES
FINAL_NCBI="08_final_NCBI/GCA_029448725.1_ASM2944872v1_genomic.fna"


# LOAD REQUIRED MODULES
module load samtools/1.15
module load bedtools/2.30.0

# 1. Seperate chromosomes from unplaced contigs

## Make a list of chromosomes 

## Unzip
#gzip -d $POLISHED 

## Index
samtools faidx $POLISHED


## Make bed and extract chr: Unplaced contigs names start with scaf
less "$POLISHED".fai | awk -F '\t' '{printf("%s\t0\t%s\n",$1,$2);}' | grep -v 'scaf' > "${POLISHED%.*}".chrs.bed


## Remove unplaced scaffolds from fasta
bedtools getfasta -fi "$POLISHED" -bed "${POLISHED%.fasta*}".chrs.bed -fullHeader > "${POLISHED%.*}".chrs.fasta

# 2. Count bases in fasta
## Count bases
CHR_BP=$(less "${POLISHED%.*}".chrs.fasta | grep -E '^[AaCcTtGgNn]+$' | perl -pe 's/[[:space:]]//g' | wc -c)
echo "$CHR_BP bases in chromosomes"


# 3. Count total bases

## Count bases
GENOME_BP=$(less $POLISHED | grep -E '^[AaCcTtGgNn]+$' | perl -pe 's/[[:space:]]//g' | wc -c)
echo "$GENOME_BP total base pairs in fasta"

# Get proportion of base pairs anchored into chromosomes
echo "$(expr $CHR_BP / $GENOME_BP) of total base pairs are anchored into chromosomes"
