#!/bin/sh

# Run on Manitou
# srun -p medium -c 10 --mem=50G --time=7-00:00:00 -J final_NCBI_RepeatClassifier -o log/final_NCBI_RepeatClassifier_%j.log /bin/sh 01_scripts/final_NCBI_RepeatClassifier.sh &


# VARIABLES
FINAL_NCBI="08_final_NCBI/GCA_029448725.1_ASM2944872v1_genomic.fna"
FINAL_NCBI_CHR="08_final_NCBI/GCA_029448725.1_ASM2944872v1_genomic_chrs.fasta"

RMOD_DIR="08_final_NCBI/RepeatModeler"
RMAS_DIR="08_final_NCBI/RepeatMasker"

STK_FILE="$RMOD_DIR/RM_117598.WedFeb71234362024/families.stk"
CONS_FILE="$RMOD_DIR/RM_117598.WedFeb71234362024/consensi.fa"

CPU=5

# LOAD REQUIRED MODULES 
module load RepeatModeler/2.0.1


# De novo detection of repeats in whole assembly
#cd "$RMOD_DIR"
#BuildDatabase -name safo ../../$FINAL_NCBI
#RepeatModeler -pa $CPU -database safo

RepeatClassifier -consensi $CONS_FILE -stockholm $STK_FILE