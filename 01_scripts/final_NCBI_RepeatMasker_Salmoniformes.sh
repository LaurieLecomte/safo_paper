#!/bin/sh

# Run on Manitou
# srun -p medium -c 10 --mem=50G --time=7-00:00:00 -J final_NCBI_RepeatMasker_Salmoniformes -o log/final_NCBI_RepeatMasker_Salmoniformes_%j.log /bin/sh 01_scripts/final_NCBI_RepeatMasker_Salmoniformes.sh &


# VARIABLES

FINAL_NCBI="08_final_NCBI/GCA_029448725.1_ASM2944872v1_genomic.fna"
FINAL_NCBI_CHR="08_final_NCBI/GCA_029448725.1_ASM2944872v1_genomic_chrs.fasta"

RMOD_DIR="08_final_NCBI/RepeatModeler"
RMAS_DIR="08_final_NCBI/RepeatMasker"

SAFO_CLASS_LIB="$RMOD_DIR/safo-families.fa"

CPU=10

# LOAD REQUIRED MODULES 
module load gnu-openmpi/4.1.4
module load exonerate/2.4.0
module load RepeatMasker/4.0.8
module load ncbiblast/2.6.0
module load python/2.7
module load maker/2.31.10


mkdir $RMAS_DIR/all_Salmoniformes

# 1. Combine the custom safo library outputted by RepeatModeler with Salmoniformes from RepBase
## Extract Salmoniformes repeat sequences
queryRepeatDatabase.pl -species 'Salmoniformes' > $RMAS_DIR/all_Salmoniformes/Salmoniformes_lib.fa

## Concatenate
cat "$RMOD_DIR/safo-families.fa" $RMAS_DIR/all_Salmoniformes/Salmoniformes_lib.fa > $RMAS_DIR/all_Salmoniformes/Salmoniformes_and_SaFo_lib.fa

# 2. Run RepeatMasker on all contigs, including unplaced contigs
#RepeatMasker -pa $CPU $FINAL_NCBI -dir $RMAS_DIR/all -gff -lib "$RMOD_DIR/safo-families.fa"
#RepeatMasker -pa $CPU $FINAL_NCBI -dir $RMAS_DIR/all_sp_salmo -species 'salmo'

RepeatMasker -pa $CPU $FINAL_NCBI -dir $RMAS_DIR/all_Salmoniformes -lib $RMAS_DIR/all_Salmoniformes/Salmoniformes_and_SaFo_lib.fa -gff