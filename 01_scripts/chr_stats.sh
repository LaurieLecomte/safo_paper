#!/bin/bash

# Get stats for chromosomes

# srun -c 1 -p small --mem=20G --time=1-00:00:00 -J chr_stats -o log/chr_stats_%j.log /bin/sh ./01_scripts/chr_stats.sh &

# VARIABLES
DRAFT_DIR="03_draft"
DRAFT="$DRAFT_DIR/flye_defaults.fasta.gz"
FINAL_NCBI="08_final_NCBI/GCA_029448725.1_ASM2944872v1_genomic.fna.gz"


# LOAD REQUIRED MODULES
module load samtools/1.15
module load bedtools/2.30.0

# 1. Seperate chromosomes from unplaced contigs

## Make a list of chromosomes 

## Unzip
gzip -d $FINAL_NCBI 

## Index
samtools faidx "${FINAL_NCBI%.*}"


## Make bed and extract chr: Unplaced contigs names start with JAQ
less "${FINAL_NCBI%.*}".fai | awk -F '\t' '{printf("%s\t0\t%s\n",$1,$2);}' | grep -v 'JAQ' > "${FINAL_NCBI%.fna*}".bed


## Remove unplaced scaffolds from fasta
bedtools getfasta -fi "${FINAL_NCBI%.*}" -bed "${FINAL_NCBI%.fna*}".bed -fullHeader > "${FINAL_NCBI%.fna*}"_chrs.fasta

# 2. Count bases in fasta
## Count bases
GENOME_BP=$(less "${FINAL_NCBI%.*}" | grep -E '^[AaCcTtGgNn]+$' | perl -pe 's/[[:space:]]//g' | wc -c)



## Extract 
#samtools faidx human_genome.fa chr1 chr2:1:2000 [...] > human_selected.fa


# 2. Count bases in fasta
## Remove unplaced scaffolds from fasta
bedtools getfasta -fi "${FINAL_NCBI%.*}" -bed $CHRS_BED -fullHeader > 03_genome/genome_no_unplaced.fasta

## Count bases
GENOME_BP=$(less 03_genome/genome_no_unplaced.fasta | grep -E '^[AaCcTtGgNn]+$' | perl -pe 's/[[:space:]]//g' | wc -c)

