#!/bin/sh

# Run RepeatMasker on chromosomes only for synteny analysis

# Run on Manitou
# srun -p medium -c 10 --mem=50G --time=7-00:00:00 -J 04.1_RepeatMasker_chrs_SaFo -o log/04.1_RepeatMasker_chrs_SaFo_%j.log /bin/sh 01_scripts/04.1_RepeatMasker_chrs_SaFo.sh &


# VARIABLES

FINAL_NCBI="08_final_NCBI/GCA_029448725.1_ASM2944872v1_genomic.fna"

SYN_DIR="synteny"

CPU=6


RMOD_DIR="08_final_NCBI/RepeatModeler"
RMAS_DIR="08_final_NCBI/RepeatMasker/all_famdb_salmonidae"

SAFO_CLASS_LIB="$RMOD_DIR/safo-families.fa"
SAFO_SALMO_LIB="$RMAS_DIR/combined_safo_Salmonidae_desc.fa" # we use the custom library produced by the final_NCBI_RepeatMasker_famdb.sh script

# LOAD REQUIRED MODULES 
module load gnu-openmpi/4.1.4
module load exonerate/2.4.0
#module load RepeatMasker/4.0.8
module load RepeatMasker/4.1.2
module load ncbiblast/2.6.0
module load python/2.7
#module load maker/2.31.10

module load samtools/1.15
module load bedtools/2.30.0
module load python/3.7

# Create synteny directory if it does not exist
if [[ ! -d $SYN_DIR ]]
then
  mkdir $SYN_DIR
fi


# 1. Extract chromosomes from fasta
## Index
#samtools faidx $FINAL_NCBI
## Make bed and extract chrs only: Remove unplaced contigs names (starting with 'scaf'), and remove mitochondrion (CM055727.1)
less "$FINAL_NCBI".fai | awk -F '\t' '{printf("%s\t0\t%s\n",$1,$2);}' | grep -E "^CM|HG|LR" | grep -v 'CM055727.1' > $SYN_DIR/SaFo.chrs_noMt.bed
less $SYN_DIR/SaFo.chrs_noMt.bed | wc -l
 
## Remove unplaced scaffolds from fasta
bedtools getfasta -fi $FINAL_NCBI -bed $SYN_DIR/SaFo.chrs_noMt.bed > $SYN_DIR/SaFo.chrs_noMt.fasta

# 2. Rename chromosomes
## Make correspondance file
samtools faidx $SYN_DIR/SaFo.chrs_noMt.fasta

if [[ -f $SYN_DIR/SaFo.chrs_noMt.corrsp.txt ]]
then
  rm $SYN_DIR/SaFo.chrs_noMt.corrsp.txt
fi


count=0
less $SYN_DIR/SaFo.chrs_noMt.fasta.fai | while read line
do
    ((count+=1))
    CHR=$(echo "$line" | cut -f1 )
    echo -e "$CHR\tChr$count" >> $SYN_DIR/SaFo.chrs_noMt.corrsp.txt 
done

## I edited line 108 of Eric's script to avoid sorting chr by new character and preserve original order. I did not change sorting command at line 117 because I know I do not have other contigs than my main chrs. 
## # rename file to prevent issues with nucmer afterwards
python3 01_scripts/rename_scaffolds.py $SYN_DIR/SaFo.chrs_noMt.fasta $SYN_DIR/SaFo.chrs_noMt.corrsp.txt 10000 $SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta


# 3. Run RepeatMasker renamed fasta
RepeatMasker -pa $CPU $SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta -dir $SYN_DIR -gff -lib $SAFO_SALMO_LIB
