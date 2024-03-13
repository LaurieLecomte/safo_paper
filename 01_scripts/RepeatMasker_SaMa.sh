#!/bin/sh

# Run on Manitou
# srun -p medium -c 6 --mem=50G --time=7-00:00:00 -J RepeatMasker_SaMa -o log/RepeatMasker_SaMa_%j.log /bin/sh 01_scripts/RepeatMasker_SaMa.sh &


# VARIABLES

FINAL_NCBI="08_final_NCBI/GCA_029448725.1_ASM2944872v1_genomic.fna"
FINAL_NCBI_CHR="08_final_NCBI/GCA_029448725.1_ASM2944872v1_genomic_chrs.fasta"

SYN_DIR="synteny/SaFo_SaMa"

SAMA="species_comparison/GCA_002910315.2_ASM291031v2_genomic.fna"
SAMA_CHRS="species_comparison/GCA_002910315.2_ASM291031v2_genomic.chrs.fasta"

RMOD_DIR="08_final_NCBI/RepeatModeler"

RMAS_DIR="08_final_NCBI/RepeatMasker"

CPU=6

SAFO_CLASS_LIB="$RMOD_DIR/safo-families.fa"
SAFO_SALMO_LIB="$RMAS_DIR/combined_safo_Salmonidae_desc.fa" # we use the custom library produced by the final_NCBI_RepeatMasker_famdb.sh scripts

# LOAD REQUIRED MODULES 
module load gnu-openmpi/4.1.4
module load exonerate/2.4.0
#module load RepeatMasker/4.0.8
module load RepeatMasker/4.1.2
module load ncbiblast/2.6.0
module load python/2.7
#module load maker/2.31.10

module load samtools/1.15
module load bedtools/2.30.0
module load python/3.7


if [[ ! -d $SYN_DIR ]]
then
  mkdir $SYN_DIR
fi

# 1. Extract chromosomes from fasta
## Index
#samtools faidx $SAMA
## Make bed and extract chr: Unplaced contigs names start with scaf, and remove mitochondrion for synteny analysis
less "$SAMA".fai | awk -F '\t' '{printf("%s\t0\t%s\n",$1,$2);}' | grep -E "^CM|HG|LR|LG" | grep -v 'CM028292.1' > $SYN_DIR/SaMa.chrs_noMt.bed
less $SYN_DIR/SaMa.chrs_noMt.bed | wc -l 
## Remove unplaced scaffolds and mitochondrion from fasta
bedtools getfasta -fi "$SAMA" -bed $SYN_DIR/SaNa.chrs_noMt.bed > $SYN_DIR/SaNa.chrs_noMt.fasta

# 2. Rename chromosomes
## Make correspondance file
samtools faidx $SYN_DIR/SaMa.chrs_noMt.fasta

if [[ -f $SYN_DIR/SaMa.chrs_noMt.corrsp.txt ]]
then
  rm $SYN_DIR/SaMa.chrs_noMt.corrsp.txt
fi


count=0
less $SYN_DIR/SaMa.chrs_noMt.fasta.fai | while read line
do
    ((count+=1))
    CHR=$(echo "$line" | cut -f1 )
    echo -e "$CHR\tChr$count" >> $SYN_DIR/SaMa.chrs_noMt.corrsp.txt 
done

## I edited line 108 of Eric's script to avoid sorting chr by new character and preserve original order. I did not change sorting command at line 117 because I know I do not have other contigs than my main chrs. 
## # rename file to prevent issues with nucmer afterwards
python3 01_scripts/rename_scaffolds.py $SYN_DIR/SaMa.chrs_noMt.fasta $SYN_DIR/SaMa.chrs_noMt.corrsp.txt 10000 $SYN_DIR/query_SaMa.chrs_noMt_renamed.fasta 

# 3. Create repeat library for Salmoniformes
queryRepeatDatabase.pl -species 'Salmoniformes' > $SYN_DIR/Salmoniformes.fa
 
# 4. Run RepeatMasker renamed fasta
RepeatMasker -pa $CPU $SYN_DIR/query_SaMa.chrs_noMt_renamed.fasta -dir $SYN_DIR -gff -lib $SAFO_SALMO_LIB
