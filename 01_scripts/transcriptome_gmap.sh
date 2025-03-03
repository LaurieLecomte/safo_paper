#!/bin/bash

# srun -p small -c 5 --mem=50G -J transcriptome_gmap -o log/transcriptome_gmap_%j.log /bin/sh 01_scripts/transcriptome_gmap.sh "species_comparison/sasaICSASG.fna" &

# VARIABLES
ASM=$1
#OUTDIR="species_comparison/transcriptome"
PROT="species_comparison/transcriptome/GCF_905237065.1_Ssal_v3.1_protein.faa.gz"
TRA="species_comparison/transcriptome/GCF_905237065.1_Ssal_v3.1_rna.fna.gz"
#INDEXED_GENOME_FOLDER="species_comparison/transcriptome/indexed_genome"

CPU=5

# LOAD REQUIRED MODULES
module load gmap/2019-03-15

# Index genome
gmap_build -d indexed_genome -D species_comparison/transcriptome/"$(basename -s '.fna' $ASM)" $ASM


zless $TRA | gmap -t $CPU \
        --dir species_comparison/transcriptome/"$(basename -s '.fna' $ASM)" \
        -d indexed_genome \
        -f gff3_gene \
        --gff3-add-separators=0 > species_comparison/transcriptome/"$(basename -s '.fna' $ASM)"/Ssal_v3.1_on_"$(basename -s '.fna' $ASM)"_gmap.gff3

    
        