#!/bin/sh

# Run on Manitou
# srun -p medium -c 10 --mem=50G --time=7-00:00:00 -J final_NCBI_RepeatModeler -o log/final_NCBI_RepeatModeler_%j.log /bin/sh 01_scripts/final_NCBI_RepeatModeler.sh &


# VARIABLES

FINAL_NCBI="08_final_NCBI/GCA_029448725.1_ASM2944872v1_genomic.fna"
FINAL_NCBI_CHR="08_final_NCBI/GCA_029448725.1_ASM2944872v1_genomic_chrs.fasta"

RM_DIR="08_final_NCBI/RepeatMasker"

SALMO="salmo/genome.fasta"

CPU=10

# LOAD REQUIRED MODULES 
module load RepeatModeler/2.0.1


# run on salmo genome
BuildDatabase -name ssalar $SALMO
RepeatModeler -pa $CPU -species 'salmo salar' $SALMO -dir salmo -gff