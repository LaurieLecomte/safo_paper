#!/bin/bash

# Get proportion of base pairs anchored into chromosomes

# srun -c 1 -p small --mem=20G --time=1-00:00:00 -J 01_stats_chr_SaFo -o log/01_stats_chr_SaFo_%j.log /bin/sh ./01_scripts/01_stats_chr_SaFo.sh &

# VARIABLES
DRAFT_DIR="03_draft"
DRAFT="$DRAFT_DIR/flye_defaults.fasta"
FINAL_NCBI="08_final_NCBI/GCA_029448725.1_ASM2944872v1_genomic.fna"


# LOAD REQUIRED MODULES
module load samtools/1.15
module load bedtools/2.30.0

# 1. Seperate chromosomes from unplaced contigs
## Index
samtools faidx $FINAL_NCBI

## Make bed and extract chr: Unplaced contigs names start with JAQ
less "$FINAL_NCBI".fai | awk -F '\t' '{printf("%s\t0\t%s\n",$1,$2);}' | grep -v 'JAQ' > "${FINAL_NCBI%.fna*}".bed

## Remove unplaced scaffolds from fasta
bedtools getfasta -fi "${FINAL_NCBI%.*}" -bed "${FINAL_NCBI%.fna*}".bed -fullHeader > "${FINAL_NCBI%.fna*}"_chrs.fasta



# 2.Get total base pairs in the final assembly
## Count bases in complete fasta
TOTAL_BP=$(less "${FINAL_NCBI%.*}" | grep -E '^[AaCcTtGgNn]+$' | perl -pe 's/[[:space:]]//g' | wc -c)
echo "$TOTAL_BP in the assembly"

# 3. Get proportion of base pairs anchored into main chromosomes
## Count bases in chr-only fasta
bedtools getfasta -fi "${FINAL_NCBI%.*}" -bed $CHRS_BED -fullHeader > 03_genome/genome_no_unplaced.fasta

## Count bases
CHR_BP=$(less 03_genome/genome_no_unplaced.fasta | grep -E '^[AaCcTtGgNn]+$' | perl -pe 's/[[:space:]]//g' | wc -c)
echo "$CHR_BP in chromosomes"
