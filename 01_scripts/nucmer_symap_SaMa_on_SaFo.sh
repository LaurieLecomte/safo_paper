#!/bin/bash

# Prepare files required for running SyMAP to identify syntenic regions
# Input is the masked reference genome , with chr names corrected to suit SyMAP (e.g. no special characters, ., _, -, ..)

# I use the conda env mummer with version 4.0.0 installed

# srun -p large -c 8 --mem=150G --time=21-00:00:00 -J nucmer_symap_SaMa_on_SaFo -o log/nucmer_symap_SaMa_on_SaFo_%j.log /bin/sh 01_scripts/nucmer_symap_SaMa_on_SaFo.sh SaMa &
# srun -p large -c 8 --mem=150G --time=21-00:00:00 -J nucmer_symap_SaSa_on_SaFo -o log/nucmer_symap_SaSa_on_SaFo_%j.log /bin/sh 01_scripts/nucmer_symap_SaSa_on_SaFo.sh SaSa &
# srun -p large -c 8 --mem=150G --time=21-00:00:00 -J nucmer_symap_SaNa_on_SaFo -o log/nucmer_symap_SaNa_on_SaFo_%j.log /bin/sh 01_scripts/nucmer_symap_SaNa_on_SaFo.sh SaNa &


# VARIABLES
CPU=8
SPECIES=$1
SYN_DIR="synteny"
OUTPUT_DIR="$SYN_DIR/SaFo_$SPECIES"

#QUERY="synteny/SaFo_"$SPECIES"/query_"$SPECIES".chrs_noMt_renamed.fasta.masked"
QUERY=$(find synten$SYN_DIR/SaFo_"$SPECIES" -type f -iname "query_$SPECIES.chrs_*renamed.fasta.masked")
REF="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta"


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

nucmer -l 100 -c 500 -t $CPU $REF $QUERY -p $OUTPUT_DIR/"$SPECIES"_on_SaFo.all

show-coords -dlTH $OUTPUT_DIR/"$SPECIES"_on_SaFo.all.delta > $OUTPUT_DIR/"$SPECIES"_on_SaFo.all.mum

