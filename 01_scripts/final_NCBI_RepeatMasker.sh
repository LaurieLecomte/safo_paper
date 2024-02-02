#!/bin/sh

# Run on Manitou
# srun -p medium -c 10 --mem=50G --time=7-00:00:00 -J final_NCBI_RepeatMasker -o log/final_NCBI_RepeatMasker_%j.log /bin/sh 01_scripts/final_NCBI_RepeatMasker.sh &


# VARIABLES

FINAL_NCBI="08_final_NCBI/GCA_029448725.1_ASM2944872v1_genomic.fna"
FINAL_NCBI_CHR="08_final_NCBI/GCA_029448725.1_ASM2944872v1_genomic_chrs.fasta"

RM_DIR="08_final_NCBI/RepeatMasker"


CPU=10

# LOAD REQUIRED MODULES 
module load exonerate/2.4.0
module load RepeatMasker/4.0.8
module load ncbiblast/2.6.0
module load python/2.7
module load gnu-openmpi/4.1.4
module load maker/2.31.10


# 1. Run RepeatMasker on primary chromosomes
RepeatMasker -pa $CPU -species 'Salvelinus fontinalis' $FINAL_NCBI_CHR -dir $RM_DIR/chr -gff -html

# 2. Run RepeatMasker on all contigs, including unplaced contigs
RepeatMasker -pa $CPU -species 'Salvelinus fontinalis' $FINAL_NCBI -dir $RM_DIR/all -gff -html