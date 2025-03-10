#!/bin/bash

# Prepare files required for running SyMAP to identify syntenic regions
# Input is the masked reference genome , with chr names corrected to suit SyMAP (e.g. no special characters, ., _, -, ..)

# I use the conda env mummer with version 4.0.0 installed

# srun -p medium -c 6 --mem=100G --time=7-00:00:00 -J test_ntsynt_unmasked_SaMa -o log/test_ntsynt_SaMa_on_SaFo_unmasked_%j.log /bin/sh 01_scripts/test_ntsynt_unmasked.sh SaMa &
# srun -p medium -c 6 --mem=100G --time=7-00:00:00 -J test_ntsynt_unmasked_SaSa -o log/test_ntsynt_SaSa_on_SaFo_unmasked_%j.log /bin/sh 01_scripts/test_ntsynt_unmasked.sh SaSa &
# srun -p medium -c 6 --mem=100G --time=7-00:00:00 -J test_ntsynt_unmasked_SaNa -o log/test_ntsynt_SaNa_on_SaFo_unmasked_%j.log /bin/sh 01_scripts/test_ntsynt_unmasked.sh SaNa &
# srun -p medium -c 6 --mem=100G --time=7-00:00:00 -J test_ntsynt_unmasked_SaFo -o log/test_ntsynt_SaFo_on_SaFo_unmasked_%j.log /bin/sh 01_scripts/test_ntsynt_unmasked.sh SaFo &


# VARIABLES
CPU=6
SPECIES=$1
SYN_DIR="synteny"
OUTPUT_DIR="$SYN_DIR/SaFo_$SPECIES"

#QUERY="synteny/SaFo_"$SPECIES"/query_"$SPECIES".chrs_noMt_renamed_unmasked.fasta.gz"


QUERY=$(find $SYN_DIR/SaFo_"$SPECIES" -type f -iname "*$SPECIES.chrs_*renamed_unmasked.fasta.gz")
#REF="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta"
#REF="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta.masked"
REF="synteny/ref_SaFo.chrs_noMt_renamed_unmasked.fasta.gz"


NTSYNT="/project/lbernatchez/users/lalec31/.conda/pkgs/ntsynt-1.0.0-py310h0dbaff4_1/bin/ntSynt"

# LOAD REQUIRED_MODULES
#module load mummer/3.23

if [[ ! -d $OUTPUT_DIR ]]
then
  mkdir $OUTPUT_DIR
fi

echo "aligning $QUERY on $REF"



$NTSYNT -d 5 -p $OUTPUT_DIR/"$SPECIES"_on_SaFo_ntSynt_unmasked $REF $QUERY

