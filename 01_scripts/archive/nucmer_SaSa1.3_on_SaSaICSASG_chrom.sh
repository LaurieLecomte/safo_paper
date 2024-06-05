#!/bin/bash

# map one query chr against all other chromosomes at once 

# This script is parallelized on each ref chr

# First extract chr in fasta for using as query 
#module load samtools/1.15
#SYN_DIR="synteny"
#OUTPUT_DIR="$SYN_DIR/SaSaICSASG_SaSa1.3"
#REF="$OUTPUT_DIR/ref_SaSaICSASG.chrs_renamed.fasta.masked"
#samtools faidx $REF
#less "$REF".fai | cut -f1 > $OUTPUT_DIR/ref_SaSaICSASG_chr_list.txt

#less $OUTPUT_DIR/ref_SaSaICSASG_chr_list.txt | while read CHR; do samtools faidx $REF $CHR > $OUTPUT_DIR/chr/r"$CHR"_ICSASG.masked.fasta; done
##samtools faidx $QUERY $CHR > $OUTPUT_DIR/q"$CHR"_renamed_masked.fasta


# parallel -a "synteny/SaSaICSASG_SaSa1.3/ref_SaSaICSASG_chr_list.txt" -j 10 srun -p small -c 2 -J r{}_SaSaICSASG_SaSa1.3_by_chrom -o log/nucmer_SaSaICSASG_on_maskedSaFo_chrom_{}_%j.log /bin/sh 01_scripts/nucmer_SaSa1.3_on_SaSaICSASG_chrom.sh {} &

# After all chr are done:

#OUTPUT_DIR="synteny/SaSaICSASG_SaSa1.3"
# file_list=$(for i in $(seq 1 29); do echo "$OUTPUT_DIR/chr_on_all/qSaSa1.3_on_rSaSaICSASG_Chr"$i".maxmatch.mum"; done)
#cat $file_list > $OUTPUT_DIR/qSaSa1.3_chrs_on_rSaSaICSASG.maxmatch.mum


# VARIABLES
CPU=2

CHR=$1
SYN_DIR="synteny"
OUTPUT_DIR="$SYN_DIR/SaSaICSASG_SaSa1.3"

REF="$OUTPUT_DIR/ref_SaSaICSASG.chrs_renamed.fasta.masked"
QUERY="$OUTPUT_DIR/query_SaSa.chrs_renamed.fasta.masked"



# LOAD REQUIRED_MODULES
#module load mummer/3.23
module load samtools/1.15


if [[ ! -d $OUTPUT_DIR ]]
then
  mkdir $OUTPUT_DIR
fi
 

if [[ ! -d $OUTPUT_DIR/chr_on_all ]]
then
  mkdir $OUTPUT_DIR/chr_on_all
fi




echo "aligning $CHR on $REF"
#using the parameters from Weissensteiner et al
# -l is the minimum length of the single match (default =20)
#-c is the minimum length of a cluster of match (defaut 65)

nucmer -t $CPU $OUTPUT_DIR/chr/r"$CHR"_ICSASG.masked.fasta $QUERY --maxmatch -p $OUTPUT_DIR/chr_on_all/qSaSa1.3_on_rSaSaICSASG_"$CHR".maxmatch

show-coords -dlTH $OUTPUT_DIR/chr_on_all/qSaSa1.3_on_rSaSaICSASG_"$CHR".maxmatch.delta > $OUTPUT_DIR/chr_on_all/qSaSa1.3_on_rSaSaICSASG_"$CHR".maxmatch.mum


