#!/bin/sh

# Run on Manitou
# srun -p large -c 15 --mem=50G --time=21-00:00:00 -J 04.3_RepeatMasker_chrs_SaSaICSASG -o log/04.3_RepeatMasker_chrs_SaSaICSASG_%j.log /bin/sh 01_scripts/04.3_RepeatMasker_chrs_SaSaICSASG.sh &


# VARIABLES

SYN_DIR="synteny/SaFo_SaSa"

SASA="species_comparison/GCA_000233375.4_ICSASG_v2_genomic.fna"

RMOD_DIR="08_final_NCBI/RepeatModeler"
RMAS_DIR="08_final_NCBI/RepeatMasker"

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


if [[ ! -d $SYN_DIR ]]
then
  mkdir $SYN_DIR
fi

# 1. Extract chromosomes from fasta
## Index
#samtools faidx $SASA
## Make bed and extract chr: Remove unplaced contigs (starting with scaf), and remove mitochondrion (no mitochondrion in Ssal genome)
less "$SASA".fai | awk -F '\t' '{printf("%s\t0\t%s\n",$1,$2);}' | grep -E "^CM|HG|LR" > $SYN_DIR/SaSa.chrs.bed
less $SYN_DIR/SaSa.chrs.bed | wc -l 
## Remove unplaced scaffolds
bedtools getfasta -fi $SASA -bed $SYN_DIR/SaSa.chrs.bed > $SYN_DIR/SaSa.chrs.fasta

# 2. Rename chromosomes
## Make correspondance file
samtools faidx $SYN_DIR/SaSa.chrs.fasta

if [[ -f $SYN_DIR/SaSa.chrs.corrsp.txt ]]
then
  rm $SYN_DIR/SaSa.chrs.corrsp.txt
fi


count=0
less $SYN_DIR/SaSa.chrs.fasta.fai | while read line
do
    ((count+=1))
    CHR=$(echo "$line" | cut -f1 )
    echo -e "$CHR\tChr$count" >> $SYN_DIR/SaSa.chrs.corrsp.txt 
done

## I edited line 108 of Eric's script to avoid sorting chr by new character and preserve original order. I did not change sorting command at line 117 because I know I do not have other contigs than my main chrs. 
## # rename file to prevent issues with nucmer afterwards
python3 01_scripts/rename_scaffolds.py $SYN_DIR/SaSa.chrs.fasta $SYN_DIR/SaSa.chrs.corrsp.txt 10000 $SYN_DIR/query_SaSa.chrs_renamed.fasta


# 3. Run RepeatMasker renamed fasta
RepeatMasker -pa $CPU $SYN_DIR/query_SaSa.chrs_renamed.fasta -dir $SYN_DIR -gff -lib $SAFO_SALMO_LIB
