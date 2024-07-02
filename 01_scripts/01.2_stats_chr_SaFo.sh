#!/bin/bash

# Get proportion of base pairs anchored into chromosomes

# srun -c 1 -p small --mem=20G --time=1-00:00:00 -J 01_stats_chr_SaFo -o log/01_stats_chr_SaFo_%j.log /bin/sh ./01_scripts/01_stats_chr_SaFo.sh &

# VARIABLES
DRAFT="03_draft/flye_defaults.fasta"
SCAF1="04_1st_scaff/BT_scaffolds_i3_FINAL.fasta"
SCAF2="05_2nd_scaf/BT_scaffolds_i3_FINAL.chromonomer.fasta"
POLISHED="06_polished/BT_scaffolds_i3_FINAL.chromonomer.hapo3x.fasta"
FINAL_NCBI="08_final_NCBI/GCA_029448725.1_ASM2944872v1_genomic.fna"


# LOAD REQUIRED MODULES
module load samtools/1.15
module load bedtools/2.30.0


for i in $SCAF2 $POLISHED $FINAL_NCBI;
do

  # 1. Seperate chromosomes from unplaced contigs  
  ## Index
  samtools faidx $i

  ## Make bed for chromosomes only: Unplaced contigs names start with scaf or JAQ
  less "$i".fai | awk -F '\t' '{printf("%s\t0\t%s\n",$1,$2);}' | grep -E -v 'scaf|JAQ' > "${i%.*}".bed

  ## Remove unplaced scaffolds from fasta
  bedtools getfasta -fi $i -bed "${i%.*}".bed -fullHeader > "${i%.*}"_chrs.fasta
  
  # 2. Count base pairs in chr-only fasta
  CHR_BP=$(less "${i%.*}"_chrs.fasta | grep -E '^[AaCcTtGgNn]+$' | perl -pe 's/[[:space:]]//g' | wc -c)
  echo "$CHR_BP in chromosomes for $i"


  # 3. Get total base pairs in the final assembly
  ## Count bases in complete fasta
  TOTAL_BP=$(less $i | grep -E '^[AaCcTtGgNn]+$' | perl -pe 's/[[:space:]]//g' | wc -c)
  echo $i
  echo "$TOTAL_BP in the assembly $i"
done


# 1. Seperate chromosomes from unplaced contigs  
  ## Index
  #samtools faidx $FINAL_NCBI

## Make bed and extract chr: Unplaced contigs names start with JAQ
#less "$FINAL_NCBI".fai | awk -F '\t' '{printf("%s\t0\t%s\n",$1,$2);}' | grep -v 'JAQ' > "${FINAL_NCBI%.fna*}".bed

## Remove unplaced scaffolds from fasta
#bedtools getfasta -fi "${FINAL_NCBI%.*}" -bed "${FINAL_NCBI%.fna*}".bed -fullHeader > "${FINAL_NCBI%.fna*}"_chrs.fasta



# 2.Get total base pairs in the final assembly
## Count bases in complete fasta
#TOTAL_BP=$(less "${FINAL_NCBI%.*}" | grep -E '^[AaCcTtGgNn]+$' | perl -pe 's/[[:space:]]//g' | wc -c)
#echo "$TOTAL_BP in the assembly"

# 3. Get proportion of base pairs anchored into main chromosomes
## Count bases in chr-only fasta
#bedtools getfasta -fi "${FINAL_NCBI%.*}" -bed $CHRS_BED -fullHeader > 03_genome/genome_no_unplaced.fasta

## Count bases
#CHR_BP=$(less 03_genome/genome_no_unplaced.fasta | grep -E '^[AaCcTtGgNn]+$' | perl -pe 's/[[:space:]]//g' | wc -c)
#echo "$CHR_BP in chromosomes"
