#!/bin/bash


# First split fasta by chr
#module load samtools/1.15
#less 02_infos/SaFo_corrected_chrs.txt | while read CHR; do samtools faidx "synteny/ref_SaFo.chrs_noMt_renamed.fasta.masked" $CHR > synteny/SaFo_"$CHR"_masked.fasta; done

# Prepare files required for running SyMAP to do self synteny
# Input is the masked reference genome , with chr names corrected to suit SyMAP (e.g. no special characters, ., _, -, ..)

# I use the conda env mummer with version 4.0.0 installed

# less 02_infos/SaFo_corrected_chrs.txt | tail -n+4 | parallel -j 20 srun -p medium -c 5 --mem=80G --time=7-00:00:00 -J {}_nucmer_self_align_chrom_vs_chrom -o log/nucmer_self_align_chrom_vs_chrom_{}_%j.log /bin/sh 01_scripts/nucmer_self_align_chrom_vs_chrom.sh {} &


# VARIABLES
CPU=5

R_CHR=$1 

SPECIES="SaFo"
SYN_DIR="synteny"
OUTPUT_DIR="$SYN_DIR/SaFo_$SPECIES"

#REF="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta"
#QUERY="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta"

GENOME="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta.masked"
#QUERY="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta.masked"

CHR_LIST="02_infos/SaFo_corrected_chrs.txt"

# LOAD REQUIRED_MODULES
#module load mummer/3.23
module load samtools/1.15


if [[ ! -d $OUTPUT_DIR ]]
then
  mkdir $OUTPUT_DIR
fi

if [[ ! -d $OUTPUT_DIR/chr_on_chr/$R_CHR ]]
then
  mkdir $OUTPUT_DIR/chr_on_chr/$R_CHR
fi


# First split fasta by chr
#module load samtools/1.15
#less $CHR_LIST | while read CHR; do samtools faidx $GENOME $CHR > $OUTPUT_DIR/chr/SaFo_"$CHR"_masked.fasta; done


# Map the reference chr R_CHR to each chr in the genome
less $CHR_LIST | while read "Q_CHR";
do
  echo "align query "$Q_CHR" on reference "$R_CHR"";
  nucmer -t $CPU $OUTPUT_DIR/chr/SaFo_"$R_CHR"_masked.fasta $OUTPUT_DIR/chr/SaFo_"$Q_CHR"_masked.fasta --maxmatch -p $OUTPUT_DIR/chr_on_chr/$R_CHR/q"$Q_CHR"_on_r"$R_CHR".maxmatch;
  show-coords -dlTH $OUTPUT_DIR/chr_on_chr/$R_CHR/q"$Q_CHR"_on_r"$R_CHR".maxmatch.delta > $OUTPUT_DIR/chr_on_chr/$R_CHR/q"$Q_CHR"_on_r"$R_CHR".maxmatch.mum;
done


# Combine outputs together
cat $OUTPUT_DIR/chr_on_chr/$R_CHR/*_on_r"$R_CHR".maxmatch.mum > $OUTPUT_DIR/chr_on_chr/$R_CHR/all_on_r"$R_CHR".maxmatch.mum

#echo "aligning $CHR on $REF"
#using the parameters from Weissensteiner et al
# -l is the minimum length of the single match (default =20)
#-c is the minimum length of a cluster of match (defaut 65)

#nucmer -t $CPU $REF $OUTPUT_DIR/q"$CHR"_renamed_masked.fasta --maxmatch --genome -p $OUTPUT_DIR/q"$CHR"_on_SaFo.maxmatch

#show-coords -dlTH $OUTPUT_DIR/q"$CHR"_on_SaFo.maxmatch.delta > $OUTPUT_DIR/q"$CHR"_on_SaFo.maxmatch.mum


