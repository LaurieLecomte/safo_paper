#!/bin/bash

# Perform self-vs-self alignents
# Nucmer crashes if all chromosomes are mapped againt all chromosomes, so I do chromosome-vs-all alignments
  
# After running aligments, I concatenate outputs together:
# After all chr are done:
#SPECIES='SaFo'
#SYN_DIR="synteny"
#OUTPUT_DIR="$SYN_DIR/SaFo_$SPECIES"
# file_list=$(for i in $(seq 1 42); do echo "synteny/SaFo_"$SPECIES"/chr_on_all/q"$SPECIES"_on_rSaFoChr"$i".maxmatch.mum"; done)
#cat $file_list > synteny/SaFo_"$SPECIES"/q"$SPECIES"chrs_on_rSaFo.maxmatch.mum


# use the conda env mummer with version 4.0.0 installed
# less "synteny/ref_SaFo.chrs_noMt_renamed.fasta.masked.fai" | cut -f1 > 02_infos/SaFo_corrected_chrs.txt 
# parallel -a 02_infos/SaFo_corrected_chrs.txt -j 20 srun -p medium -c 8 --mem=100G --time=7-00:00:00 -J {}_nucmer_self_chr_vs_all -o log/nucmer_symap_self_align_chrom_{}_%j.log /bin/sh 01_scripts/nucmer_symap_self_align_chrom.sh {} &


# VARIABLES
CPU=8

CHR=$1
SPECIES="SaFo"
SYN_DIR="synteny"
OUTPUT_DIR="$SYN_DIR/SaFo_$SPECIES"

REF="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta.masked"

# LOAD REQUIRED_MODULES
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
less $CHR_LIST | while read CHR; do samtools faidx $GENOME $CHR > $OUTPUT_DIR/chr/SaFo_"$CHR"_masked.fasta; done
samtools faidx $QUERY $CHR > $OUTPUT_DIR/q"$CHR"_renamed_masked.fasta

# Map each chr against all chromosomes
echo "aligning $CHR on $REF"
#using the parameters from Weissensteiner et al
# -l is the minimum length of the single match (default =20)
#-c is the minimum length of a cluster of match (defaut 65)

nucmer -t $CPU $REF $OUTPUT_DIR/chr/SaFo_"$CHR"_masked.fasta --maxmatch -p $OUTPUT_DIR/chr_on_all/q"$CHR"_on_SaFo.maxmatch

show-coords -dlTH $OUTPUT_DIR/chr_on_all/q"$CHR"_on_SaFo.maxmatch.delta > $OUTPUT_DIR/chr_on_all/q"$CHR"_on_SaFo.maxmatch.mum


