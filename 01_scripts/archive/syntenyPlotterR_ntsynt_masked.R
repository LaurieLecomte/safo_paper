SaNa_synblocks <- read.delim("~/SaFo_paper/SaNa_on_SaFo_masked_ntSynt.synteny_blocks.tsv", header=FALSE,
                            col.names = c('block_num', 'file', 'chr', 'start', 'end', 'strand', 'minimizers', 'disc'))
SaNa_synblocks$type <- sapply(X = SaNa_synblocks$file, FUN = function(x) unlist(strsplit(x, split = '_'))[1])

SaNa_synblocks <- select(SaNa_synblocks, ! 'file')

library(dplyr)

pivot_wider(SaNa_synblocks,
            names_from = 'block_num',
            values_from = 'type')

SaNa_synblocks <- 
pivot_wider(SaNa_synblocks,
            names_from = 'type',
            values_from = c('chr', 'start', 'end', 'strand'))
SaNa_synblocks$size_ref <- SaNa_synblocks$end_ref - SaNa_synblocks$start_ref
SaNa_synblocks$size_query <- SaNa_synblocks$end_query - SaNa_synblocks$start_query


#SaNa_synblocks <- subset(SaNa_synblocks, size_ref > 2000000 & size_query > 2000000)


SaNa_blocks_fmt <- SaNa_synblocks[, c('chr_ref', 'start_ref', 'end_ref', 'chr_query', 'start_query', 'end_query')]

#SaNa_blocks_fmt$orient <- '+'
SaNa_blocks_fmt$orient <- SaNa_synblocks$strand_ref

SaNa_blocks_fmt$rID <- 'SaFo'
SaNa_blocks_fmt$qID <- 'SaNa'
SaNa_blocks_fmt$chr_ref <- sapply(X = SaNa_blocks_fmt$chr_ref, FUN = function(x) substr(x, start = 4, stop = 6))
SaNa_blocks_fmt$chr_query <- sapply(X = SaNa_blocks_fmt$chr_query, FUN = function(x) substr(x, start = 4, stop = 6))



write.table(SaNa_blocks_fmt, file = "~/SaFo_paper/ntsynt_SaNa_SaFo_masked_blocks_fmt.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")



SaNa_chr_order <- c(4,3,2,1,8,5,7,16,21,26,6,9,12,10,15,28,11,20,19,29,31,27,13,23,33,25,24,18,17,35,14,32,34,41,30,36,22,37,38,39,40,42)
SaNa_lengths <- q1chr_lengths[order(match(q1chr_lengths$chrNum, SaNa_chr_order)),]



SaSa_synblocks <- read.delim("~/SaFo_paper/SaSa_on_SaFo_masked_ntSynt.synteny_blocks.tsv", header=FALSE,
                             col.names = c('block_num', 'file', 'chr', 'start', 'end', 'strand', 'minimizers', 'disc'))
SaSa_synblocks$type <- sapply(X = SaSa_synblocks$file, FUN = function(x) unlist(strsplit(x, split = '_'))[1])

SaSa_synblocks <- select(SaSa_synblocks, ! 'file')

library(dplyr)

pivot_wider(SaSa_synblocks,
            names_from = 'block_num',
            values_from = 'type')

SaSa_synblocks <- 
  pivot_wider(SaSa_synblocks,
              names_from = 'type',
              values_from = c('chr', 'start', 'end', 'strand'))
SaSa_synblocks$size_ref <- SaSa_synblocks$end_ref - SaSa_synblocks$start_ref
SaSa_synblocks$size_query <- SaSa_synblocks$end_query - SaSa_synblocks$start_query


#SaSa_synblocks <- subset(SaSa_synblocks, size_ref > 2000000 & size_query > 2000000)


SaSa_blocks_fmt <- SaSa_synblocks[, c('chr_ref', 'start_ref', 'end_ref', 'chr_query', 'start_query', 'end_query')]

SaSa_blocks_fmt$orient <- SaSa_synblocks$strand_ref
SaSa_blocks_fmt$rID <- 'SaFo'
SaSa_blocks_fmt$qID <- 'SaSa'
SaSa_blocks_fmt$chr_ref <- sapply(X = SaSa_blocks_fmt$chr_ref, FUN = function(x) substr(x, start = 4, stop = 6))
SaSa_blocks_fmt$chr_query <- sapply(X = SaSa_blocks_fmt$chr_query, FUN = function(x) substr(x, start = 4, stop = 6))


write.table(SaSa_blocks_fmt, file = "~/SaFo_paper/ntsynt_SaSa_SaFo_masked_blocks_fmt.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

SaSa_chr_order <- c(1,13,12,24,23,5,16,29,13,10,11,7,16,
                  4,10,1,
                  9, 14, 22, 21,15,1,
                  27,25,20,19,3,6,20,9,28,15,14,9,6,11,18,18,2,17)

SaSa_lengths <- q2chr_lengths[order(match(q2chr_lengths$chrNum, unique(SaSa_chr_order))),]


# 4. Test syntenyPlotteR --------------------------------------------------
all_chrlengths <- rbind(SaNa_lengths, rchr_lengths[, 1:3], SaSa_lengths)

write.table(all_chrlengths, file = "~/SaFo_paper/all_chrlengths.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")



library(syntenyPlotteR)
draw_linear_laurie(directory = "~/SaFo_paper",
            paste0('test_ntsynt_masked', format(Sys.time(), "%Y%m%d.%H%M%S")),
            sizefile = "~/SaFo_paper/all_chrlengths.txt", 
            "~/SaFo_paper/ntsynt_SaNa_SaFo_masked_blocks_fmt.txt",
            "~/SaFo_paper/ntsynt_SaSa_SaFo_masked_blocks_fmt.txt",
            colours = rep(c('black', 'blue'), (nrow(all_chrlengths)+1)/2))



draw_linear_laurie(
  directory = "~/SaFo_paper",
  paste0('test', format(Sys.time(), "%Y%m%d.%H%M%S")),
  sizefile = "~/SaFo_paper/all_chrlengths.txt", 
  "~/SaFo_paper/querySaNa_to_refSaFo/results/blocks_fmt.txt",
  "~/SaFo_paper/querySaSa_to_refSaFo/results/blocks_fmt.txt",
  colours = rep(c('grey20', 'blue'), (nrow(all_chrlengths)+1)/2))



#q2chr_lengths <- q2chr_lengths[order(match(q2chr_lengths$chr, unique(q2_blocks_fmt$qChrCM))),]
q2chr_lengths <- q2chr_lengths[order(match(q2chr_lengths$chr, unique(q2_blocks_fmt$qChr))),]




