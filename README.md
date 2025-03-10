# Post-assembly analysis of the brook trout (*Salvelinus fontinalis*) genome assembly ASM2944872v1

This repository contains the scripts used for various post-assembly analyses described in the paper Chromosome-level genome assembly of a doubled haploid brook trout (*Salvelinus fontinalis*) (in production). 

The folder `01_scripts` contains scripts used for post-assembly steps (stats, assembly comparison, ...). These scripts generate files that are meant to be outputted in the `species_comparison` and `synteny` destination folders.


The ASM2944872v1 assembly was produced from a combination of ccs mode PacBio long reads (120X), high coverage Illumina NovaSeq6000 short reads (70-80X) and Hi-C data, and is available from GenBank (accession GCA_029448725.1) and RefSeq (accession [GCF_029448725.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_029448725.1/)).

For details on assembly procedure, tools and parameters, see the Methods section in the paper.

These analyses relied mainly on the final assembly fasta downloaded from RefSeq, but also on the intermediary files produced during the assembly procedure and the following genome assemblies of other salmonid species:
* Lake trout (*Salvelinus namaycush*): SaNama_1.0; [GCF_016432855.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_016432855.1/) (SaNa)
* Dolly Varden (*Salvelinus* sp. IW2-2015): ASM291031v2; [GCF_002910315.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002910315.2/) (referred to as SaMa as *S. malma*)
* Rainbow trout (*Oncorhynchus mykiss*): USDA_OmykA_1.1; [GCF_013265735.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_013265735.2/) (OMyk)
* River trout (*Salmo trutta*): fSalTru1.1; [GCF_901001165.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_901001165.1/) (SaTr)
* Atlantic salmon (*Salmo salar*): ICSASG_v2; [GCF_000233375.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000233375.1/) (SaSa)
* Lake whitefish (*Coregonus clupeaformis*): ASM2061545v1; [GCF_020615455.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_020615455.1/) (CoCl)

For synteny analyses, only the brook trout (SaFo), lake trout (SaNa) and Atlantic salmon (SaSa) genomes were used.


## Overview 
* Compute various assembly statistics and metrics for the whole SaFo assembly (chromosomes + contigs): `01.1_assembly_stats_SaFo.sh`
* Compute statistics and metrics for **main chromosomes** only: `01.2_stats_chr_SaFo.sh`
* Characterize repeats in the whole SaFo assembly
    * De novo repeat detection and generation of custom repeat library: `02.1_final_RepeatModeler.sh`
    * Repeat masking using the custom and the Salmonidae library: `02.2_final_RepeatMasker.sh`
* Compute statistics for various salmonid genome assemblies for comparison purposes: `03.0_species_comp.sh`
* Explore synteny with related species SaNa and SaSa
    * Repeat masking **chromosomes** of SaFo, SaNa and SaSa using the custom and the Salmonidae library: `04.1_RepeatMasker_chrs_SaFo.sh`, `04.2_RepeatMasker_chrs_SaNa.sh` and `04.3_RepeatMasker_chrs_SaSaICSASG.sh`
    * Homologous block detection between SaFo-SaNa and SaFo-SaSa: `04.4_homology_ntsynt.sh`
* Explore self synteny and homeologous regions within SaFo chromosomes
    * Self vs self genome alignment of SaFo chromosomes: `05.1_nucmer_self_align_chrom.sh`
    * Estimate % homology between homeologous regions: `05.2_lastz_align_blocks.sh`
    * Estimate coverage across the genome to highligh putatively collapsed regions: `05.3_depth.sh`
    
Shield: [![CC BY-SA 4.0][cc-by-sa-shield]][cc-by-sa]

This work is licensed under a
[Creative Commons Attribution-ShareAlike 4.0 International License][cc-by-sa].

[![CC BY-SA 4.0][cc-by-sa-image]][cc-by-sa]

[cc-by-sa]: http://creativecommons.org/licenses/by-sa/4.0/
[cc-by-sa-image]: https://licensebuttons.net/l/by-sa/4.0/88x31.png
[cc-by-sa-shield]: https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg