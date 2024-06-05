#!/bin/bash

# I use the conda env sibeliaz

# srun -p medium -c 6 --mem=100G --time=7-00:00:00 -J test_sibeliaz_SaMa -o log/test_sibeliaz_SaMa_on_SaFo_%j.log /bin/sh 01_scripts/test_sibeliaz.sh SaMa &
# srun -p medium -c 6 --mem=100G --time=7-00:00:00 -J test_sibeliaz_SaSa -o log/test_sibeliaz_SaSa_on_SaFo_%j.log /bin/sh 01_scripts/test_sibeliaz.sh SaSa &
# srun -p medium -c 6 --mem=100G --time=7-00:00:00 -J test_sibeliaz_SaNa -o log/test_sibeliaz_SaNa_on_SaFo_%j.log /bin/sh 01_scripts/test_sibeliaz.sh SaNa &
# srun -p medium -c 6 --mem=100G --time=7-00:00:00 -J test_sibeliaz_SaFo -o log/test_sibeliaz_SaFo_on_SaFo_%j.log /bin/sh 01_scripts/test_sibeliaz.sh SaFo &


# VARIABLES
CPU=6
SPECIES=$1
SYN_DIR="synteny"
OUTPUT_DIR="$SYN_DIR/SaFo_$SPECIES"

#QUERY="synteny/SaFo_SaSa/query_SaSa.chrs_renamed.masked.fasta"
QUERY=$(find $SYN_DIR/SaFo_"$SPECIES" -type f -iname "query_$SPECIES.chrs_*renamed.fasta.masked")
#REF="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta"
REF="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta.masked"

# LOAD REQUIRED_MODULES
#module load mummer/3.23

if [[ ! -d $OUTPUT_DIR ]]
then
  mkdir $OUTPUT_DIR
fi

echo "aligning $QUERY on $REF"


#sibeliaz -t $CPU -o $OUTPUT_DIR/"$SPECIES"_on_SaFo_masked_sibeliaz $REF $QUERY

maf2synteny $OUTPUT_DIR/"$SPECIES"_on_SaFo_masked_sibeliaz/alignment.maf -o $OUTPUT_DIR/"$SPECIES"_on_SaFo_masked_sibeliaz