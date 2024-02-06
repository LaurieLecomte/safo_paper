#/bin/bash

# Test assembly-stats 

# srun -p small -c 1 -J test_assembly-stats -o log/test_assembly-stats_%j.log /bin/sh 01_scripts/test_assembly_stats.sh &

# VARIABLES
#DRAFT_DIR="03_draft"
DRAFT="03_draft/flye_defaults.fasta"
SCAF1="04_1st_scaff/BT_scaffolds_i3_FINAL.fasta"
SCAF2="05_2nd_scaf/BT_scaffolds_i3_FINAL.chromonomer.fasta"
POLISHED="06_polished/BT_scaffolds_i3_FINAL.chromonomer.hapo3x.fasta"
FINAL_NCBI="08_final_NCBI/GCA_029448725.1_ASM2944872v1_genomic.fna"

# LOAD REQUIRED MODULES

#gzip -d $DRAFT_DIR/flye_defaults.fasta.gz

#assembly-stats -t $DRAFT_DIR/flye_defaults.fasta
 
#assembly-stats -t $DRAFT $SCAF1 $SCAF2 $POLISHED $FINAL_NCBI


# test BBmap
module load BBTools

stats.sh $DRAFT > "${DRAFT%.*}".stats.txt

stats.sh $SCAF1 > "${SCAF1%.*}".stats.txt

stats.sh $SCAF2 > "${SCAF2%.*}".stats.txt

stats.sh $POLISHED > "${POLISHED%.*}".stats.txt

stats.sh $FINAL_NCBI > "${FINAL_NCBI%.*}".stats.txt
