#/bin/bash

# Compute basic assembly stats for different the different SaFo assembly versions using the assembly-stats utility and the stats.sh script from BBTools
# First load a conda env where assembly stats was installed

# srun -p small -c 1 -J 01_assembly_stats_SaFo -o log/01_assembly_stats_SaFo_%j.log /bin/sh 01_scripts/01_assembly_stats_SaFo.sh &

# VARIABLES
#DRAFT_DIR="03_draft"
DRAFT="03_draft/flye_defaults.fasta"
SCAF1="04_1st_scaff/BT_scaffolds_i3_FINAL.fasta"
SCAF2="05_2nd_scaf/BT_scaffolds_i3_FINAL.chromonomer.fasta"
POLISHED="06_polished/BT_scaffolds_i3_FINAL.chromonomer.hapo3x.fasta"
FINAL_NCBI="08_final_NCBI/GCA_029448725.1_ASM2944872v1_genomic.fna"

# LOAD REQUIRED MODULES
module load BBTools/36.92

# Run assembly-stats on each version
assembly-stats -t $DRAFT $SCAF1 $SCAF2 $POLISHED $FINAL_NCBI


# Run stats.sh (from BBTools) on each version
stats.sh $DRAFT > "${DRAFT%.*}".BBTools.stats.txt

stats.sh $SCAF1 > "${SCAF1%.*}".BBTools.stats.txt

stats.sh $SCAF2 > "${SCAF2%.*}".BBTools.stats.txt

stats.sh $POLISHED > "${POLISHED%.*}".BBTools.stats.txt

stats.sh $FINAL_NCBI > "${FINAL_NCBI%.*}".BBTools.stats.txt
