#!/bin/bash

# Prepare files required for running SyMAP to identify syntenic regions, because it cannot handle self-alignements for huge genomes. We first map the genome to itself with nucmer, then we trick symap to use the .mum file outputted by nucmer.
# Input is the masked reference genome outputted by RepeatMasker (script utils/RepeatMasker_manitou.sh), with chr names corrected to suit SyMAP (e.g. no special characters, ., _, -, ..)
# Trials resulted in nucmer failure due to a problematic region on chr OV354449, so we had to "patch" this region by masking it with Ns prior to running nucmer.

# We first "patch" the problematic region in masked reference genome, then run nucmer. The .mum file was then transfered on the machine on which SyMAP is installed. 

# srun -p medium -c 8 --mem=100G --time=7-00:00:00 -J mummer_self_synteny -o log/nucmer_self_synteny_%j.log /bin/sh 01_scripts/nucmer_symap_maxmatch_self.sh &


# VARIABLES
CPU=6
#SPECIES=$1
SYN_DIR="synteny"
OUTPUT_DIR="$SYN_DIR"

#QUERY="synteny/SaFo_"$SPECIES"/query_"$SPECIES".chrs_noMt_renamed.fasta.masked"
#QUERY=$(find $SYN_DIR/SaFo_"$SPECIES" -type f -iname "query_$SPECIES.chrs_*renamed.fasta.masked")
QUERY="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta"
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

nucmer -t $CPU $REF $QUERY --maxmatch --genome -p $OUTPUT_DIR/SaFo_on_SaFo.maxmatch

show-coords -dlTH $OUTPUT_DIR/SaFo_on_SaFo.maxmatch.delta > $OUTPUT_DIR/SaFo_on_SaFo.maxmatch.mum

