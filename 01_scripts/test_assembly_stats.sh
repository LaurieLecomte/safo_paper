#/bin/bash

# Test assembly-stats 

# srun -p small -c 1 -J test_assembly-stats -o log/test_assembly-stats_%j.log /bin/sh 01_scripts/test_assembly_stats.sh &

# VARIABLES
DRAFT_DIR="03_draft"
#DRAFT="03_draft/flye_defaults.fasta.gz"
#SCAF1
#SCAF2
#SCAF3
#POLISHED
#FINAL=
#FINAL_NCBI

# LOAD REQUIRED MODULES

gzip -d $DRAFT_DIR/flye_defaults.fasta.gz

assembly-stats -t $DRAFT_DIR/flye_defaults.fasta
 
 