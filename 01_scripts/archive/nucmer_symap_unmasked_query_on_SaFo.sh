#!/bin/bash

# Prepare files required for running SyMAP to identify syntenic regions
# Input is the masked reference genome , with chr names corrected to suit SyMAP (e.g. no special characters, ., _, -, ..)

# I use the conda env mummer with version 4.0.0 installed

# srun -p medium -c 6 --mem=100G --time=7-00:00:00 -J nucmer_symap_unmasked_SaMa_on_SaFo -o log/nucmer_symap_unmasked_SaMa_on_SaFo_%j.log /bin/sh 01_scripts/nucmer_symap_unmasked_query_on_SaFo.sh SaMa &
# srun -p medium -c 6 --mem=100G --time=7-00:00:00 -J nucmer_symap_unmasked_SaSa_on_SaFo -o log/nucmer_symap_unmasked_SaSa_on_SaFo_%j.log /bin/sh 01_scripts/nucmer_symap_unmasked_query_on_SaFo.sh SaSa &
# srun -p medium -c 6 --mem=100G --time=7-00:00:00 -J nucmer_symap_unmasked_SaNa_on_SaFo -o log/nucmer_symap_unmasked_SaNa_on_SaFo_%j.log /bin/sh 01_scripts/nucmer_symap_unmasked_query_on_SaFo.sh SaNa &


# VARIABLES
CPU=6
SPECIES=$1
SYN_DIR="synteny"
OUTPUT_DIR="$SYN_DIR/SaFo_$SPECIES"

#QUERY="synteny/SaFo_"$SPECIES"/query_"$SPECIES".chrs_noMt_renamed.fasta.masked"
QUERY=$(find $SYN_DIR/SaFo_"$SPECIES" -type f -iname "query_$SPECIES.chrs_*renamed.fasta")
REF="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta"
#REF="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta.masked"

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

nucmer -l 100 -c 500 --genome --maxmatch -t $CPU $REF $QUERY -p $OUTPUT_DIR/"$SPECIES"_on_SaFo.maxmatch.unmasked

show-coords -dlTH $OUTPUT_DIR/"$SPECIES"_on_SaFo.maxmatch.unmasked.delta > $OUTPUT_DIR/"$SPECIES"_on_SaFo.maxmatch.unmasked.all.mum

