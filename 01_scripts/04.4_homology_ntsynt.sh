#!/bin/bash

# Run ntsynt on species pair to find homology between chromosomes
# Run in a conda env where ntsynt is installed

# parallel -j 2 srun -p small -c 4 -J 04.4_homology_ntsynt_{} -o log/04.4_homology_ntsynt_{}_%j.log /bin/sh 01_scripts/04.4_homology_ntsynt.sh {} & ::: SaSaICSASG SaNa



# VARIABLES
CPU=4
SPECIES=$1
SYN_DIR="synteny"
OUTPUT_DIR="$SYN_DIR/SaFo_$SPECIES"

QUERY=$(find $SYN_DIR/SaFo_"$SPECIES" -type f -iname "query_$SPECIES.chrs_*renamed.fasta.masked")
#REF="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta" 
REF="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta.masked"

NTSYNT="/project/lbernatchez/users/lalec31/.conda/pkgs/ntsynt-1.0.0-py310h0dbaff4_1/bin/ntSynt"

# LOAD REQUIRED_MODULES
#module load mummer/3.23

if [[ ! -d $OUTPUT_DIR ]]
then
  mkdir $OUTPUT_DIR
fi

echo "running ntsynt with $QUERY on $REF"


#module load mash/2.0

#mash dist $QUERY $REF

$NTSYNT -d 10 --block_size 1000000 --indel 100000 --merge 1000000 --w_rounds 500 250 -t $CPU -p $OUTPUT_DIR/"$SPECIES"_on_SaFo_ntSynt_b1M_d10 $REF $QUERY
#$NTSYNT -d 1 --block_size 100000 --indel 10000 --merge 100000 --w_rounds 100 10 -t $CPU -p $OUTPUT_DIR/"$SPECIES"_on_SaFo_ntSynt_b1M_m1M_d1 $REF $QUERY
