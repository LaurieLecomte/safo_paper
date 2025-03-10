#!/bin/sh

# Run on Manitou
# srun -p medium -c 10 --mem=50G --time=7-00:00:00 -J final_NCBI_RepeatMasker_famdb -o log/final_NCBI_RepeatMasker_famdb_%j.log /bin/sh 01_scripts/final_NCBI_RepeatMasker_famdb.sh &


# VARIABLES

FINAL_NCBI="08_final_NCBI/GCA_029448725.1_ASM2944872v1_genomic.fna"
FINAL_NCBI_CHR="08_final_NCBI/GCA_029448725.1_ASM2944872v1_genomic_chrs.fasta"

RMOD_DIR="08_final_NCBI/RepeatModeler"
RMAS_DIR="08_final_NCBI/RepeatMasker/all_famdb_salmonidae"

#SAFO_CLASS_LIB="08_final_NCBI/RepeatModeler/RM_117598.WedFeb71234362024/consensi.fa.classified"
SAFO_CLASS_LIB="$RMOD_DIR/safo-families.fa"

CPU=10

# LOAD REQUIRED MODULES 
module load gnu-openmpi/4.1.4
module load exonerate/2.4.0
module load RepeatMasker/4.1.2
module load ncbiblast/2.6.0
module load python/2.7
#module load maker/2.31.10


if [[ -f "$RMAS_DIR/safo-families_renamed.fa" ]]
then
  rm "$RMAS_DIR/safo-families_renamed.fa"
fi

# 1. Create custom lib using the famdb utility
famdb.py  families --format 'fasta_name' --descendants 'Salmonidae' --include-class-in-name > $RMAS_DIR/famdb_Salmonidae_desc.fa

# 2. Combine with custom lib produced by 
## First edit sequence header in safo lib
less $SAFO_CLASS_LIB | while read line; do echo $line | sed -E 's/(^.+)\ \(.+/\1\ \@Salvelinus_fontinalis/' >> "$RMAS_DIR/safo-families_renamed.fa" ; done

## Concatenate
cat $RMAS_DIR/famdb_Salmonidae_desc.fa "$RMAS_DIR/safo-families_renamed.fa" > $RMAS_DIR/combined_safo_Salmonidae_desc.fa

# 3. Run RepeatMasker on all contigs, using the combined custom library
RepeatMasker -pa $CPU $FINAL_NCBI -dir $RMAS_DIR -gff -lib $RMAS_DIR/combined_safo_Salmonidae_desc.fa
