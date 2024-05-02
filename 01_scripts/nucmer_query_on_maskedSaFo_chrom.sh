#!/bin/bash

# map one query chr against all other chromosomes at once 

# Prepare files required for running SyMAP to do self synteny
# Input is the masked reference genome , with chr names corrected to suit SyMAP (e.g. no special characters, ., _, -, ..)

# I use the conda env mummer with version 4.0.0 installed
# -a 02_infos/SaFo_corrected_chrs.txt
# parallel -a 02_infos/SaFo_corrected_chrs.txt -j 20 srun -p medium -c 8 --mem=100G --time=7-00:00:00 -J {}_nucmer_query_on_maskedSaFo_chrom -o log/nucmer_query_on_maskedSaFo_chrom_{}_%j.log /bin/sh 01_scripts/nucmer_query_on_maskedSaFo_chrom.sh {} SaNa &


# VARIABLES
CPU=2

CHR=$1
SPECIES=$2
SYN_DIR="synteny"
OUTPUT_DIR="$SYN_DIR/SaFo_$SPECIES"

#REF="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta"
#QUERY="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta"


QUERY=$(find $SYN_DIR/SaFo_"$SPECIES" -type f -iname "query_$SPECIES.chrs_*renamed.fasta.masked")
REF="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta.masked"
#QUERY="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta.masked"

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


# Extract chr in fasta for using as query 
# First split fasta by chr
#module load samtools/1.15
#less $CHR_LIST | while read CHR; do samtools faidx $GENOME $CHR > $OUTPUT_DIR/chr/SaFo_"$CHR"_masked.fasta; done
##samtools faidx $QUERY $CHR > $OUTPUT_DIR/q"$CHR"_renamed_masked.fasta



echo "aligning $CHR on $REF"
#using the parameters from Weissensteiner et al
# -l is the minimum length of the single match (default =20)
#-c is the minimum length of a cluster of match (defaut 65)

nucmer -t $CPU $SYN_DIR/SaFo_SaFo/chr/SaFo_"$CHR"_masked.fasta $QUERY --maxmatch -p $OUTPUT_DIR/chr_on_all/q"$SPECIES"_on_rSaFo"$CHR".maxmatch

show-coords -dlTH $OUTPUT_DIR/chr_on_all/q"$SPECIES"_on_rSaFo"$CHR".maxmatch.delta > $OUTPUT_DIR/chr_on_all/q"$SPECIES"_on_rSaFo"$CHR".maxmatch.mum


