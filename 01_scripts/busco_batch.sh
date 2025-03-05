#!/bin/bash

# srun -p medium --time=0-12:00:00 -c 20 --mem=100G -J busco_batch -o log/busco_batch_%j.log /bin/sh 01_scripts/busco_batch.sh &

# VARIABLES
BUSCO_DIR="species_comparison/busco"

#SAFO="$BUSCO_DIR/in/GCF_029448725.1_ASM2944872v1_genomic.fna"
#SANA="$BUSCO_DIR/GCF_016432855.1_SaNama_1.0_genomic.fna"
#OMYK="$BUSCO_DIR/GCF_013265735.2_USDA_OmykA_1.1_genomic.fna"
#SATR="$BUSCO_DIR/GCF_901001165.1_fSalTru1.1_genomic.fna"
#SASA="$BUSCO_DIR/GCF_905237065.1_Ssal_v3.1_genomic.fna"
#COCL="$BUSCO_DIR/GCF_020615455.1_ASM2061545v1_genomic.fna"
#SAMA="$BUSCO_DIR/GCF_002910315.2_ASM291031v2_genomic.fna"

#SASA_IC="$BUSCO_DIR/GCF_000233375.1_ICSASG_v2_genomic.fna"   #ICSASG assembly v2

CPU=20

# LOAD REQUIRED MODULES
module load busco/5.8.2
source activate busco-5.8.2

busco --in $BUSCO_DIR --lineage_dataset actinopterygii_odb10 --mode genome --cpu $CPU --out_path $BUSCO_DIR -f
#busco --in $BUSCO_DIR/inputs/GCF_029448725.1_ASM2944872v1_gnomon_protein.faa --lineage_dataset actinopterygii_odb10 --mode protein --cpu $CPU --out_path $BUSCO_DIR