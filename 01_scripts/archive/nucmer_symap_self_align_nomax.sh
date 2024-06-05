#!/bin/bash

# Prepare files required for running SyMAP to do self synteny
# Input is the masked reference genome , with chr names corrected to suit SyMAP (e.g. no special characters, ., _, -, ..)

# I use the conda env mummer with version 4.0.0 installed

# srun -p medium -c 12 --mem=100G --time=7-00:00:00 -J nucmer_symap_self_align_nomax -o log/nucmer_symap_self_align_nomax_%j.log /bin/sh 01_scripts/nucmer_symap_self_align_nomax.sh SaFo &


# VARIABLES
CPU=12
SPECIES=$1
SYN_DIR="synteny"
OUTPUT_DIR="$SYN_DIR/SaFo_$SPECIES"

#REF="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta"
#QUERY="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta"

REF="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta.masked"
QUERY="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta.masked"


# LOAD REQUIRED_MODULES
#module load mummer/3.23

if [[ ! -d $OUTPUT_DIR ]]
then
  mkdir $OUTPUT_DIR
fi


echo "aligning $QUERY on $REF"
#using the parameters from Weissensteiner et al
# -l is the minimum length of the single match (default =20)
#-c is the minimum length of a cluster of match (defaut 65)

nucmer -t $CPU $REF $QUERY --genome -p $OUTPUT_DIR/"$SPECIES"_on_SaFo_nomax

show-coords -dlTH $OUTPUT_DIR/"$SPECIES"_on_SaFo_nomax.delta > $OUTPUT_DIR/"$SPECIES"_on_SaFo_nomax.mum


