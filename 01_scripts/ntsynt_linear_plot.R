# linear plot for ntsynt results
library(dplyr)
library(tidyr)

REF_SIZES <- "~/SaFo_paper/SaFo.chrs_noMt.chrsizes.txt"
REF_CHRS <- "~/SaFo_paper/SaFo.chrs_for_blocks.txt"


Q1_BLOCKS <- "~/SaFo_paper/ntsynt/SaNa_on_SaFo_ntSynt_b1M_d10.synteny_blocks.tsv"
Q1_SIZES <- "~/SaFo_paper/querySaNa_to_refSaFo/SaNa.chrs_noMt.chrsizes.txt"
Q1_CHRS <- "~/SaFo_paper/querySaNa_to_refSaFo/SaNa.chrs_for_blocks.txt"


Q2_BLOCKS <- "~/SaFo_paper/ntsynt/SaSa_on_SaFo_ntSynt_b1M_d10.synteny_blocks.tsv"
Q2_SIZES <- "~/SaFo_paper/querySaSa_to_refSaFo/SaSa.chrs.chrsizes.txt"
Q2_CHRS <- "~/SaFo_paper/querySaSa_to_refSaFo/SaSa.chrs_for_blocks.txt"

MIN_BLOCK_SIZE <- 1000000

# 1. Import files for reference SaFo --------------------------------------
rchrsizes <- read.delim(REF_SIZES, header=FALSE,
                        col.names = c('chr', 'size'))

rchrs_for_blocks <- read.delim(REF_CHRS, header=FALSE,
                               col.names = c('chrCM', 'chrSymap'))
rchrs_for_blocks$chrNum <- 1:nrow(rchrs_for_blocks)

rchr_lengths <- merge(rchrsizes, rchrs_for_blocks, by.x = 'chr', by.y = 'chrCM', sort = FALSE)
rchr_lengths <- rchr_lengths[, c('chrNum', 'size')]
rchr_lengths$ID <- 'SaFo'

# Add BC linkage groups 
LG_crrsp <- read.delim("~/SaFo_paper/linkage_group_correspondance_Claire_Ben.txt")

LG_crrsp$Ben_BC <- paste0('BC', LG_crrsp$Ben_BC)

rchr_lengths <- merge(x = rchr_lengths, y =LG_crrsp, by.x = 'chrNum', by.y = 'Claire_assembly')
#rchr_lengths$Ben_BC <- paste0('BC', rchr_lengths$Ben_BC)
rchr_lengths$chrNum <- paste0(rchr_lengths$chrNum, ' (', rchr_lengths$Ben_BC, ')')
rchr_lengths <- rchr_lengths[, c('chrNum', 'size', 'ID')]


# 2. Import files for query 1 ---------------------------------------------
q1_synblocks <- read.delim(Q1_BLOCKS, header=FALSE,
                             col.names = c('block_num', 'file', 'chr', 'start', 'end', 'strand', 'minimizers', 'disc'))
q1_synblocks$type <- sapply(X = q1_synblocks$file, FUN = function(x) unlist(strsplit(x, split = '_'))[1])

q1_synblocks <- select(q1_synblocks, ! 'file')

q1_synblocks <- 
  pivot_wider(q1_synblocks,
              names_from = 'type',
              values_from = c('chr', 'start', 'end', 'strand'))
q1_synblocks$size_ref <- q1_synblocks$end_ref - q1_synblocks$start_ref
q1_synblocks$size_query <- q1_synblocks$end_query - q1_synblocks$start_query


#q1_synblocks <- subset(q1_synblocks, size_ref > MIN_BLOCK_SIZE & size_query > MIN_BLOCK_SIZE)


q1_synblocks_fmt <- q1_synblocks[, c('chr_ref', 'start_ref', 'end_ref', 'chr_query', 'start_query', 'end_query')]

#q1_synblocks_fmt$orient <- '+'
q1_synblocks_fmt$orient <- q1_synblocks$strand_ref

q1_synblocks_fmt$rID <- 'SaFo'
q1_synblocks_fmt$qID <- 'SaNa'
#q1_synblocks_fmt$chr_ref <- sapply(X = q1_synblocks_fmt$chr_ref, FUN = function(x) paste0('Ssa', substr(x, start = 4, stop = 6)))
#q1_synblocks_fmt$chr_query <- sapply(X = q1_synblocks_fmt$chr_query, FUN = function(x) paste0('Sna', substr(x, start = 4, stop = 6)))
q1_synblocks_fmt$chr_ref <- sapply(X = q1_synblocks_fmt$chr_ref, FUN = function(x) substr(x, start = 4, stop = 6))
q1_synblocks_fmt$chr_query <- sapply(X = q1_synblocks_fmt$chr_query, FUN = function(x) substr(x, start = 4, stop = 6))


# Add BC from Sutherland et al. 2016
q1_synblocks_fmt <- merge(x = q1_synblocks_fmt, y = LG_crrsp, by.x = 'chr_ref', by.y = 'Claire_assembly')

q1_synblocks_fmt$chr_ref <- paste0(q1_synblocks_fmt$chr_ref, ' (', q1_synblocks_fmt$Ben_BC, ')')




#write.table(q1_synblocks_fmt[, 1:9], file = paste0(unlist(strsplit(Q1_BLOCKS, '.tsv')), '_BC_fmt.txt'),
write.table(q1_synblocks_fmt[, 1:9], file = paste0(unlist(strsplit(Q1_BLOCKS, '.tsv')), '_fmt.txt'),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

q1chrsizes <- read.delim(Q1_SIZES, header=FALSE,
                         col.names = c('chr', 'size'))

#less "synteny/SaFo_SaSa/SaSa.chrs.corrsp.txt" | sed -E 's/^(HG.+)\:0\-[0-9]+/\1/' > "synteny/SaFo_SaSa/SaSa.chrs_for_blocks.txt"
q1chrs_for_blocks <- read.delim(Q1_CHRS, header=FALSE,
                                col.names = c('chrHG', 'chrSymap'))

q1chrs_for_blocks$chrNum <- 1:nrow(q1chrs_for_blocks)


q1chr_lengths <- merge(q1chrsizes, q1chrs_for_blocks, by.x = 'chr', by.y = 'chrHG', sort = FALSE)
q1chr_lengths <- q1chr_lengths[, c('chrNum', 'size')]
q1chr_lengths$ID <- 'SaNa'


q1_chr_order <- c(4,3,2,1,8,5,7,16,21,26,6,9,12,10,15,28,11,20,19,29,31,27,13,23,33,25,24,18,17,35,14,32,34,41,30,36,22,37,38,39,40,42)
q1chr_lengths <- q1chr_lengths[order(match(q1chr_lengths$chrNum, q1_chr_order)),]



# 3. Import files for query 2 ---------------------------------------------
q2_synblocks <- read.delim(Q2_BLOCKS, header=FALSE,
                           col.names = c('block_num', 'file', 'chr', 'start', 'end', 'strand', 'minimizers', 'disc'))
q2_synblocks$type <- sapply(X = q2_synblocks$file, FUN = function(x) unlist(strsplit(x, split = '_'))[1])

q2_synblocks <- select(q2_synblocks, ! 'file')

q2_synblocks <- 
  pivot_wider(q2_synblocks,
              names_from = 'type',
              values_from = c('chr', 'start', 'end', 'strand'))
q2_synblocks$size_ref <- q2_synblocks$end_ref - q2_synblocks$start_ref
q2_synblocks$size_query <- q2_synblocks$end_query - q2_synblocks$start_query


#q2_synblocks <- subset(q2_synblocks, size_ref > MIN_BLOCK_SIZE & size_query > MIN_BLOCK_SIZE)


q2_synblocks_fmt <- q2_synblocks[, c('chr_ref', 'start_ref', 'end_ref', 'chr_query', 'start_query', 'end_query')]

#q2_synblocks_fmt$orient <- '+'
q2_synblocks_fmt$orient <- q2_synblocks$strand_ref

q2_synblocks_fmt$rID <- 'SaFo'
q2_synblocks_fmt$qID <- 'SaSa'
#q2_synblocks_fmt$chr_ref <- sapply(X = q2_synblocks_fmt$chr_ref, FUN = function(x) paste0('Ssa', substr(x, start = 4, stop = 6)))
#q2_synblocks_fmt$chr_query <- sapply(X = q2_synblocks_fmt$chr_query, FUN = function(x) paste0('Sna', substr(x, start = 4, stop = 6)))
q2_synblocks_fmt$chr_ref <- sapply(X = q2_synblocks_fmt$chr_ref, FUN = function(x) substr(x, start = 4, stop = 6))
q2_synblocks_fmt$chr_query <- sapply(X = q2_synblocks_fmt$chr_query, FUN = function(x) substr(x, start = 4, stop = 6))


# Add BC from Sutherland et al. 2016
q2_synblocks_fmt <- merge(x = q2_synblocks_fmt, y = LG_crrsp, by.x = 'chr_ref', by.y = 'Claire_assembly')

q2_synblocks_fmt$chr_ref <- paste0(q2_synblocks_fmt$chr_ref, ' (', q2_synblocks_fmt$Ben_BC, ')')


#write.table(q2_synblocks_fmt[, 1:9], file = paste0(unlist(strsplit(Q2_BLOCKS, '.tsv')), '_BC_fmt.txt'),
write.table(q2_synblocks_fmt[, 1:9], file = paste0(unlist(strsplit(Q2_BLOCKS, '.tsv')), '_fmt.txt'),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


q2chrsizes <- read.delim(Q2_SIZES, header=FALSE,
                         col.names = c('chr', 'size'))

#less "synteny/SaFo_SaSa/SaSa.chrs.corrsp.txt" | sed -E 's/^(HG.+)\:0\-[0-9]+/\1/' > "synteny/SaFo_SaSa/SaSa.chrs_for_blocks.txt"
q2chrs_for_blocks <- read.delim(Q2_CHRS, header=FALSE,
                                col.names = c('chrHG', 'chrSymap'))

q2chrs_for_blocks$chrNum <- 1:nrow(q2chrs_for_blocks)


q2chr_lengths <- merge(q2chrsizes, q2chrs_for_blocks, by.x = 'chr', by.y = 'chrHG', sort = FALSE)
q2chr_lengths <- q2chr_lengths[, c('chrNum', 'size')]
q2chr_lengths$ID <- 'SaSa'

#q2_chr_order <- c(4,3,2,1,8,5,7,16,21,26,6,9,12,10,15,28,11,20,19,29,31,27,13,23,33,25,24,18,17,35,14,32,34,41,30,36,22,37,38,39,40,42)
#q2chr_lengths <- q1chr_lengths[order(match(q2chr_lengths$chrNum, q2_chr_order)),]
q2_chr_order <- c(1,13,12,24,23,5,16,29,13,10,11,7,16,
                  4,10,1,
                  9, 14, 22, 21,15,1,
                  27,25,20,19,3,6,20,9,28,15,14,9,6,11,18,18,2,17)

q2chr_lengths <- q2chr_lengths[order(match(q2chr_lengths$chrNum, unique(q2_chr_order))),]

# 4. Test syntenyPlotteR --------------------------------------------------
all_chrlengths <- rbind(q1chr_lengths, rchr_lengths, q2chr_lengths)
                                     
write.table(all_chrlengths, file = "~/SaFo_paper/all_chrlengths.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


library(syntenyPlotteR)
# Source custom version of draw_linear() function
source("custom_draw_linear.R")
                                     
draw_linear_laurie(
  directory = "~/SaFo_paper",
  paste0('ntsynt_b1M_d10_', format(Sys.time(), "%Y%m%d.%H%M%S")),
  #sizefile = "~/SaFo_paper/all_chrlengths_BC.txt", 
  sizefile = "~/SaFo_paper/all_chrlengths.txt", 
  paste0(unlist(strsplit(Q1_BLOCKS, '.tsv')), '_fmt.txt'),
  paste0(unlist(strsplit(Q2_BLOCKS, '.tsv')), '_fmt.txt'),
  colours = rep(c('grey20', 'blue'), nrow(rchr_lengths))
  )
