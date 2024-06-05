#!/bin/bash

# Test Syri using nucmer's output
# I use the conda env mummer with version 4.0.0 installed

# srun -p medium -c 6 --time=7-00:00 -J test_syri -o log/test_syri_%j.log /bin/sh 01_scripts/test_syri.sh SaMa &



# VARIABLES
CPU=6
SPECIES=$1
SYN_DIR="synteny"
OUTPUT_DIR="$SYN_DIR/SaFo_$SPECIES"

#QUERY="synteny/SaFo_"$SPECIES"/query_"$SPECIES".chrs_noMt_renamed.fasta.masked"
QUERY=$(find $SYN_DIR/SaFo_"$SPECIES" -type f -iname "query_$SPECIES.chrs_*renamed.fasta.masked")
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

#nucmer -l 100 -c 500 --genome -t $CPU $REF $QUERY -p $OUTPUT_DIR/"$SPECIES"_on_SaFo.all

#show-coords -dlTH $OUTPUT_DIR/"$SPECIES"_on_SaFo.all.delta > $OUTPUT_DIR/"$SPECIES"_vs_SaFo.all.mum



delta-filter -m -i 90 -l 100 $OUTPUT_DIR/"$SPECIES"_on_SaFo.all.delta > $OUTPUT_DIR/"$SPECIES"_on_SaFo.all.filt.delta     # Remove small and lower quality alignments
show-coords -THrd $OUTPUT_DIR/"$SPECIES"_on_SaFo.all.filt.delta > $OUTPUT_DIR/"$SPECIES"_on_SaFo.all.filt.coords          # Convert alignment information to a .TSV format as required by SyRI

syri -c $OUTPUT_DIR/"$SPECIES"_on_SaFo.all.filt.coords -d $OUTPUT_DIR/"$SPECIES"_on_SaFo.all.filt.delta -r $REF -q $QUERY --nc $CPU --prefix $OUTPUT_DIR/"$SPECIES"_on_SaFo.all.syri
#python3 $PATH_TO_PLOTSR $OUTPUT_DIR/"$SPECIES"_on_SaFo.all.syri $REF $QUERY -H 8 -W 5



#nucmer --maxmatch -c 100 -b 500 -l 50 $REF $QUERY -p $OUTPUT_DIR/"$SPECIES"_on_SaFo.maxmatch                                        # Whole genome alignment. Any other alignment can also be used.
#delta-filter -m -i 90 -l 100 $OUTPUT_DIR/"$SPECIES"_on_SaFo.maxmatch.delta > $OUTPUT_DIR/"$SPECIES"_on_SaFo.maxmatch.filt.delta     # Remove small and lower quality alignments
#show-coords -THrd $OUTPUT_DIR/"$SPECIES"_on_SaFo.maxmatch.filt.delta > $OUTPUT_DIR/"$SPECIES"_on_SaFo.maxmatch.filt.coords          # Convert alignment information to a .TSV format as required by SyRI

#python3 $PATH_TO_SYRI -c $OUTPUT_DIR/"$SPECIES"_on_SaFo.maxmatch.filt.coords -d $OUTPUT_DIR/"$SPECIES"_on_SaFo.maxmatch.filt.delta -r $REF -q $QUERY --nc $CPU --prefix $OUTPUT_DIR/"$SPECIES"_on_SaFo.maxmatch.syri
#python3 $PATH_TO_PLOTSR $OUTPUT_DIR/"$SPECIES"_on_SaFo.maxmatch.syri $REF $QUERY -H 8 -W 5