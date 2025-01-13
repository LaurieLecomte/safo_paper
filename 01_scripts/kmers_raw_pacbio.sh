#!/bin/bash

# Estimate coverage and genome size from raw long reads using kmers

# srun -p small -c 4 --mem=100G -J kmers_raw_pacbio -o log/kmers_raw_pacbio_%j.log /bin/sh 01_scripts/kmers_raw_pacbio.sh 21 &


# VARIABLES
RAW_LR="raw_data/all_ccs_fontinalis.fq.gz"


CPU=4

KMER_SIZE=$1

# LOAD REQUIRED MODULES
module load jellyfish/2.3.0


# Unzip 
#gzip -c -d $RAW_LR > "${RAW_LR%.*}"

# Count kmers
jellyfish count --mer-len $KMER_SIZE --size 160M --threads $CPU --canonical "${RAW_LR%.*}" --output "${RAW_LR%%.*}"_counts_k"$KMER_SIZE".jf 


# Get freq histogram
jellyfish histo --threads $CPU "${RAW_LR%%.*}"_counts_k"$KMER_SIZE".jf > "${RAW_LR%%.*}"_k"$KMER_SIZE".histo

