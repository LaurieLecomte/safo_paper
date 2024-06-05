#!/bin/bash

# Prepare files required for running SyMAP to identify syntenic regions
# Input is the masked reference genome , with chr names corrected to suit SyMAP (e.g. no special characters, ., _, -, ..)

# I use the conda env mummer with version 4.0.0 installed

# srun -p medium -c 6 --mem=100G --time=7-00:00:00 -J nucmer_softmask_SaMa_on_SaFo -o log/nucmer_softmask_SaMa_on_SaFo_%j.log /bin/sh 01_scripts/nucmer_softmask_query_on_SaFo.sh SaMa &
# srun -p medium -c 6 --mem=100G --time=7-00:00:00 -J nucmer_softmask_SaSa_on_SaFo -o log/nucmer_softmask_SaSa_on_SaFo_%j.log /bin/sh 01_scripts/nucmer_softmask_query_on_SaFo.sh SaSa &
# srun -p medium -c 6 --mem=100G --time=7-00:00:00 -J nucmer_softmask_SaNa_on_SaFo -o log/nucmer_softmask_SaNa_on_SaFo_%j.log /bin/sh 01_scripts/nucmer_softmask_query_on_SaFo.sh SaNa &


# VARIABLES
CPU=6
SPECIES=$1
SYN_DIR="synteny"
OUTPUT_DIR="$SYN_DIR/SaFo_$SPECIES"

#QUERY="synteny/SaFo_"$SPECIES"/query_"$SPECIES".chrs_noMt_renamed.fasta.masked"
QUERY_UM=$(find $SYN_DIR/SaFo_"$SPECIES" -type f -iname "query_$SPECIES.chrs_*renamed.fasta")
QUERY=$(find $SYN_DIR/SaFo_"$SPECIES" -type f -iname "query_$SPECIES.chrs_*renamed.fasta.masked")

REF="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta"  # softmasked ref
#REF="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta.masked"

# LOAD REQUIRED_MODULES
#module load mummer/3.23
module load seqkit/2.8.0
module load bedtools/2.31.1
 
 
if [[ ! -d $OUTPUT_DIR ]]
then
  mkdir $OUTPUT_DIR
fi


# Convert hardmasked query to softmasked
seqkit locate --bed -rPp "N+" $QUERY > $SYN_DIR/SaFo_"$SPECIES"/"$SPECIES"_hardmasked_regions.bed

bedtools maskfasta -soft -fi $QUERY_UM -bed $SYN_DIR/SaFo_"$SPECIES"/"$SPECIES"_hardmasked_regions.bed -fo $SYN_DIR/SaFo_"$SPECIES"/"$SPECIES"_renamed_softmasked.fasta

echo "aligning softmasked $SYN_DIR/SaFo_"$SPECIES"/"$SPECIES"_renamed_softmasked.fasta on softmasked $REF"
#using the parameters from Weissensteiner et al
# -l is the minimum length of the single match (default =20)
#-c is the minimum length of a cluster of match (defaut 65)

nucmer -l 100 -c 500 --genome -t $CPU $REF $SYN_DIR/SaFo_"$SPECIES"/"$SPECIES"_renamed_softmasked.fasta -p $OUTPUT_DIR/"$SPECIES"_on_SaFo.softmasked

show-coords -dlTH $OUTPUT_DIR/"$SPECIES"_on_SaFo.softmasked.delta > $OUTPUT_DIR/"$SPECIES"_on_SaFo.softmasked.mum

