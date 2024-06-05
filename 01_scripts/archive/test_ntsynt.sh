#!/bin/bash

# Prepare files required for running SyMAP to identify syntenic regions
# Input is the masked reference genome , with chr names corrected to suit SyMAP (e.g. no special characters, ., _, -, ..)

# I use the conda env mummer with version 4.0.0 installed

# srun -p medium -c 6 --mem=100G --time=7-00:00:00 -J test_ntsynt_SaMa -o log/test_ntsynt_SaMa_on_SaFo_%j.log /bin/sh 01_scripts/test_ntsynt.sh SaMa &
# srun -p medium -c 6 --mem=100G --time=7-00:00:00 -J test_ntsynt_SaSa -o log/test_ntsynt_SaSa_on_SaFo_%j.log /bin/sh 01_scripts/test_ntsynt.sh SaSa &
# srun -p medium -c 6 --mem=100G --time=7-00:00:00 -J test_ntsynt_SaNa -o log/test_ntsynt_SaNa_on_SaFo_%j.log /bin/sh 01_scripts/test_ntsynt.sh SaNa &

# srun -p small -c 4 --mem=50G --time=1-00:00:00 -J test_ntsynt_SaSaICSASG -o log/test_ntsynt_SaSaICSASG_on_SaFo_%j.log /bin/sh 01_scripts/test_ntsynt.sh SaSaICSASG &

# srun -p small -c 4 --mem=50G --time=1-00:00:00 -J test_ntsynt_OMyk -o log/test_ntsynt_OMyk_on_SaFo_%j.log /bin/sh 01_scripts/test_ntsynt.sh OMyk &

# srun -p medium -c 6 --mem=100G --time=7-00:00:00 -J test_ntsynt_SaFo -o log/test_ntsynt_SaFo_on_SaFo_%j.log /bin/sh 01_scripts/test_ntsynt.sh SaFo &


# VARIABLES
CPU=4
SPECIES=$1
SYN_DIR="synteny"
OUTPUT_DIR="$SYN_DIR/SaFo_$SPECIES"

#QUERY="synteny/SaFo_"$SPECIES"/query_"$SPECIES".chrs_noMt_renamed.fasta.masked"
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

echo "aligning $QUERY on $REF"


#module load mash/2.0

#mash dist $QUERY $REF

#$NTSYNT -d 10 --block_size 1000000 --indel 100000 --merge 1000000 --w_rounds 500 250 -t $CPU -p $OUTPUT_DIR/"$SPECIES"_on_SaFo_ntSynt_b1M_d10 $REF $QUERY
$NTSYNT -d 1 --block_size 100000 --indel 10000 --merge 100000 --w_rounds 100 10 -t $CPU -p $OUTPUT_DIR/"$SPECIES"_on_SaFo_ntSynt_b1M_m1M_d1 $REF $QUERY
