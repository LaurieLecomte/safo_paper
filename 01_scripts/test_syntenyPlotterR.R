# Test 


# 1. Import files for reference SaFo --------------------------------------
rchrsizes <- read.delim("~/SaFo_paper/SaFo.chrs_noMt.chrsizes.txt", header=FALSE,
                        col.names = c('chr', 'size'))

rchrs_for_blocks <- read.delim("~/SaFo_paper/SaFo.chrs_for_blocks.txt", header=FALSE,
                               col.names = c('chrCM', 'chrSymap'))
rchrs_for_blocks$chrNum <- 1:nrow(rchrs_for_blocks)

rchr_lengths <- merge(rchrsizes, rchrs_for_blocks, by.x = 'chr', by.y = 'chrCM', sort = FALSE)
rchr_lengths <- rchr_lengths[, c('chrNum', 'size')]
rchr_lengths$ID <- 'SaFo'





# 2. Import files for query SaNa ------------------------------------------
# Alignment
blocks <- read.delim("~/SaFo_paper/querySaNa_to_refSaFo/results/blocks",
                     col.names = c('qChr', 'rChr', 'block', 'qStart', 'qEnd', 'rStart', 'rEnd', 'hits', 'qGenes', 'rGenes', 'qPctGenes', 'rPctGenes'))


blocks_fmt <- blocks[, c('rChr', 'rStart', 'rEnd', 'qChr', 'qStart', 'qEnd')]
blocks_fmt$orient <- '+'
blocks_fmt$rID <- 'SaFo'
blocks_fmt$qID <- 'SaNa'


write.table(blocks_fmt, file = "~/SaFo_paper/querySaNa_to_refSaFo/results/blocks_fmt.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

#less "synteny/SaFo_SaNa/SaNa.chrs_noMt.fasta.fai" | cut -f1,2 | sed -E 's/^(CM.+)\:0\-[0-9]+/\1/'  > "synteny/SaFo_SaNa/SaNa.chrs_noMt.chrsizes.txt"
q1chrsizes <- read.delim("~/SaFo_paper/querySaNa_to_refSaFo/SaNa.chrs_noMt.chrsizes.txt", header=FALSE,
                       col.names = c('chr', 'size'))

#less "synteny/SaFo_SaNa/SaNa.chrs_noMt.corrsp.txt" | sed -E 's/^(CM.+)\:0\-[0-9]+/\1/' > "synteny/SaFo_SaNa/SaNa.chrs_for_blocks.txt"
q1chrs_for_blocks <- read.delim("~/SaFo_paper/querySaNa_to_refSaFo/SaNa.chrs_for_blocks.txt", header=FALSE,
                                   col.names = c('chrCM', 'chrSymap'))

q1chrs_for_blocks$chrNum <- 1:nrow(q1chrs_for_blocks)


q1chr_lengths <- merge(q1chrsizes, q1chrs_for_blocks, by.x = 'chr', by.y = 'chrCM', sort = FALSE)
q1chr_lengths <- q1chr_lengths[, c('chrNum', 'size')]
q1chr_lengths$ID <- 'SaNa'

#all_chrlengths <- rbind(q1chr_lengths, rchr_lengths)
#write.table(all_chrlengths, file = "~/SaFo_paper/querySaNa_to_refSaFo/SaNa_SaFo.chrlengths.txt",
#            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


# 3. Import files for query SaMa ------------------------------------------
# Alignment
blocks <- read.delim("~/SaFo_paper/querySaMa_to_refSaFo/results/blocks",
                     col.names = c('qChr', 'rChr', 'block', 'qStart', 'qEnd', 'rStart', 'rEnd', 'hits', 'qGenes', 'rGenes', 'qPctGenes', 'rPctGenes'))


blocks_fmt <- blocks[, c('rChr', 'rStart', 'rEnd', 'qChr', 'qStart', 'qEnd')]
blocks_fmt$orient <- '+'
blocks_fmt$rID <- 'SaFo'
blocks_fmt$qID <- 'SaMa'


write.table(blocks_fmt, file = "~/SaFo_paper/querySaMa_to_refSaFo/results/blocks_fmt.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

#less "synteny/SaFo_SaMa/SaMa.chrs_noMt.fasta.fai" | cut -f1,2 | sed -E 's/^(CM.+)\:0\-[0-9]+/\1/'  > "synteny/SaFo_SaMa/SaMa.chrs_noMt.chrsizes.txt"
q2chrsizes <- read.delim("~/SaFo_paper/querySaMa_to_refSaFo/SaMa.chrs_noMt.chrsizes.txt", header=FALSE,
                        col.names = c('chr', 'size'))

#less "synteny/SaFo_SaMa/SaMa.chrs_noMt.corrsp.txt" | sed -E 's/^(CM.+)\:0\-[0-9]+/\1/' > "synteny/SaFo_SaMa/SaMa.chrs_for_blocks.txt"
q2chrs_for_blocks <- read.delim("~/SaFo_paper/querySaMa_to_refSaFo/SaMa.chrs_for_blocks.txt", header=FALSE,
                               col.names = c('chrCM', 'chrSymap'))

q2chrs_for_blocks$chrNum <- 1:nrow(q2chrs_for_blocks)


q2chr_lengths <- merge(q2chrsizes, q2chrs_for_blocks, by.x = 'chr', by.y = 'chrCM', sort = FALSE)
q2chr_lengths <- q2chr_lengths[, c('chrNum', 'size')]
q2chr_lengths$ID <- 'SaMa'


#write.table(all_chrlengths, file = "~/SaFo_paper/querySaMa_to_refSaFo/SaMa_SaFo.chrlengths.txt",
#            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")



# 4. Test syntenyPlotteR --------------------------------------------------
all_chrlengths <- rbind(q1chr_lengths, rchr_lengths, q2chr_lengths)

write.table(all_chrlengths, file = "~/SaFo_paper/all_chrlengths.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")



library(syntenyPlotteR)
draw.linear(directory = "~/SaFo_paper",
            'test',
            sizefile = "~/SaFo_paper/all_chrlengths.txt", 
            "~/SaFo_paper/querySaNa_to_refSaFo/results/blocks_fmt.txt",
            "~/SaFo_paper/querySaMa_to_refSaFo/results/blocks_fmt.txt",
            colours = rep(c('black', 'blue'), (nrow(all_chrlengths)+1)/2))
