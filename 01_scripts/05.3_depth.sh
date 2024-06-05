#!/bin/bash

# Estimate shprt read depth per window 
# Use path to mapped short reads file as first positional argument

# srun -p small -c 4 --time=1-00:00:00 -J 05.3_depth -o log/05.3_depth_%j.log /bin/sh 01_scripts/05.3_depth.sh "wgs_sample_preparation/06_aligned/SaFo_cat_1.trimmed.sorted.bam" &


# VARIABLES
CPU=4
SYN_DIR="synteny"
OUTPUT_DIR="$SYN_DIR/SaFo_SaFo"

REF="$SYN_DIR/ref_SaFo.chrs_noMt_renamed.fasta.masked"

WIN_SIZE=1000000

BAM=$1
SAMPLE=$(basename -s '.trimmed.sorted.bam' $BAM)

SYN_BLOCKS="$OUTPUT_DIR/blocks.bed"


# LOAD REQUIRED MODULES
module load mosdepth/0.3.6

WIN_MB=$(( $WIN_SIZE / 1000000 ))



# Run modepth to get coverage by window
mosdepth -t $CPU -n --fast-mode --by $WIN_SIZE $OUTPUT_DIR/mosdepth_"$SAMPLE"_"$WIN_MB"Mb $BAM


# Get overlap between depth bed and syntenic blocks bed
zless $OUTPUT_DIR/mosdepth_"$SAMPLE"_"$WIN_MB"Mb.regions.bed.gz | sed -E 's/^Chr([0-9]+)/\1/' > $OUTPUT_DIR/mosdepth_"$SAMPLE"_"$WIN_MB"Mb.regions_chrnum.bed
bedtools intersect -a $OUTPUT_DIR/mosdepth_"$SAMPLE"_"$WIN_MB"Mb.regions_chrnum.bed -b $SYN_BLOCKS -wao > $OUTPUT_DIR/intersect_depth_synblocks_"$SAMPLE"_"$WIN_MB"Mb.txt


# Run modepth to get average coverage for whole genome
mosdepth -t $CPU -n --fast-mode --by $WIN_SIZE $OUTPUT_DIR/mosdepth_"$SAMPLE"_"$WIN_MB"Mb_whole $BAM