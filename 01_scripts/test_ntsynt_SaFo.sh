#!/bin/bash

# srun -p medium -c 6 --mem=100G --time=7-00:00:00 -J test_ntsynt_SaFo -o log/test_ntsynt_SaFo_on_SaFo_%j.log /bin/sh 01_scripts/test_ntsynt_SaFo.sh &


# VARIABLES
CPU=6
SPECIES='SaFo'
SYN_DIR="synteny"
OUTPUT_DIR="$SYN_DIR/SaFo_$SPECIES"

QUERY="synteny/SaFo_"$SPECIES"/query_"$SPECIES".chrs_noMt_renamed.fasta.masked"
#QUERY=$(find $SYN_DIR/SaFo_"$SPECIES" -type f -iname "query_$SPECIES.chrs_*renamed.fasta.masked")
#REF="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta"
REF="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta.masked"

NTSYNT="/project/lbernatchez/users/lalec31/.conda/pkgs/ntsynt-1.0.0-py310h0dbaff4_1/bin/ntSynt"

# LOAD REQUIRED_MODULES
#module load mummer/3.23

if [[ ! -d $OUTPUT_DIR ]]
then
  mkdir $OUTPUT_DIR
fi

echo "aligning $QUERY on $REF"



$NTSYNT -d 5 --block_size 500 --indel 10000 --merge 10000 --w_rounds 100 10 -t $CPU -p $OUTPUT_DIR/"$SPECIES"_on_SaFo_d5_ntSynt $REF $QUERY

#--block_size 1000 --indel 50000 --merge 100000 --w_rounds 250 100 
#--block_size 500 --indel 10000 --merge 10000 --w_rounds 100 10