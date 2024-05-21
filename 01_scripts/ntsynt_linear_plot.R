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
#all_chrlengths <- rbind(rchr_lengths, q2chr_lengths)

#write.table(all_chrlengths, file = "~/SaFo_paper/all_chrlengths_BC.txt",
write.table(all_chrlengths, file = "~/SaFo_paper/all_chrlengths.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

#write.table(all_chrlengths, file = "~/SaFo_paper/Q2_REF_chrlengths.txt",
 #           row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")



library(syntenyPlotteR)

draw_linear_laurie(
  directory = "~/SaFo_paper",
  paste0('ntsynt_b1M_d10_', format(Sys.time(), "%Y%m%d.%H%M%S")),
  #sizefile = "~/SaFo_paper/all_chrlengths_BC.txt", 
  sizefile = "~/SaFo_paper/all_chrlengths.txt", 
  #paste0(unlist(strsplit(Q1_BLOCKS, '.tsv')), '_BC_fmt.txt'),
  #paste0(unlist(strsplit(Q2_BLOCKS, '.tsv')), '_BC_fmt.txt'),
  paste0(unlist(strsplit(Q1_BLOCKS, '.tsv')), '_fmt.txt'),
  paste0(unlist(strsplit(Q2_BLOCKS, '.tsv')), '_fmt.txt'),
  colours = rep(c('grey20', 'blue'), nrow(rchr_lengths))
  #colours = rep(c('grey20', 'blue'), (nrow(all_chrlengths)+1)/2)
  )





chr_blocks <- (unique(rchr_lengths$chrNum))
cols_blocks <- vector(mode = 'character', length = nrow(rchr_lengths))
#hex_cols <- (viridisLite::viridis(n = length(chr_blocks), option = 'D'))
#hex_cols <- c(viridisLite::viridis(n = length(chr_blocks)/2, option = 'D'), viridisLite::viridis(n = length(chr_blocks)/2, option = 'C'))
col_pair <- rep(c('grey20', 'blue'), (nrow(rchr_lengths)+1)/2)

for (i in 1:length(chr_blocks)) {
  names(cols_blocks)[i] <- chr_blocks[i]
  cols_blocks[i] <- col_pair[i]
}

cols_vec_Q1 <- cols_blocks[match(q1_synblocks_fmt$chr_ref, names(cols_blocks))]
cols_vec_Q2 <- cols_blocks[match(q2_synblocks_fmt$chr_ref, names(cols_blocks))]

full_cols_vec <- c(cols_vec_Q1, cols_blocks, cols_vec_Q2)

draw_linear_laurie(
  directory = "~/SaFo_paper",
  paste0('ntsynt_b1M_d10_', format(Sys.time(), "%Y%m%d.%H%M%S")),
  #sizefile = "~/SaFo_paper/all_chrlengths_BC.txt", 
  sizefile = "~/SaFo_paper/all_chrlengths.txt", 
  #paste0(unlist(strsplit(Q1_BLOCKS, '.tsv')), '_BC_fmt.txt'),
  #paste0(unlist(strsplit(Q2_BLOCKS, '.tsv')), '_BC_fmt.txt'),
  paste0(unlist(strsplit(Q1_BLOCKS, '.tsv')), '_fmt.txt'),
  paste0(unlist(strsplit(Q2_BLOCKS, '.tsv')), '_fmt.txt'),
  colours = full_cols_vec)
  #colours = rep(c('grey20', 'blue'), (nrow(all_chrlengths)+1)/2)
)




q2_synblocks_dist <- 
  q2_synblocks %>% 
  group_by(chr_ref) %>%  
  arrange(chr_query, chr_ref, start_ref) %>% 
  mutate(distance = (start_ref - lag(end_ref, n = 1, default = 0)))

q2_synblocks_dist$orient <- q2_synblocks_dist$strand_ref

q2_synblocks_dist$rID <- 'SaFo'
q2_synblocks_dist$qID <- 'SaSa'

q2_synblocks_dist$chr_ref <- sapply(X = q2_synblocks_dist$chr_ref, FUN = function(x) substr(x, start = 4, stop = 6))
q2_synblocks_dist$chr_query <- sapply(X = q2_synblocks_dist$chr_query, FUN = function(x) substr(x, start = 4, stop = 6))


# Add BC from Sutherland et al. 2016
#q2_synblocks_dist_fmt <- merge(x = q2_synblocks_dist, y = LG_crrsp, by.x = 'chr_ref', by.y = 'Claire_assembly')
#q2_synblocks_dist_fmt$chr_ref <- paste0(q2_synblocks_dist_fmt$chr_ref, ' (', q2_synblocks_dist_fmt$Ben_BC, ')')

q2_synblocks_dist_rfmt <- q2_synblocks_dist

for (i in 1:nrow(q2_synblocks_dist)){
  if (q2_synblocks_dist$distance[i] < 100000){
    
    if (q2_synblocks_dist$orient[i] == q2_synblocks_dist$orient[i-1] & q2_synblocks_dist$chr_ref[i] == q2_synblocks_dist$chr_ref[i-1] & q2_synblocks_dist$chr_query[i] == q2_synblocks_dist$chr_query[i-1]){
      print('pouf')
      q2_synblocks_dist_rfmt$end_ref[i-1] <- q2_synblocks_dist$end_ref[i]
      q2_synblocks_dist_rfmt$end_query[i-1] <- q2_synblocks_dist$end_query[i]
    }
    
  }
}

q2_synblocks_dist_rfmt <- subset(q2_synblocks_dist_rfmt, distance >= 100000)
#q2_synblocks_dist_rfmt <- subset(q2_synblocks_dist_rfmt, size_ref > 2000000)

q2_synblocks_dist_rfmt <- subset(q2_synblocks_dist_rfmt, size_ref > MIN_BLOCK_SIZE & size_query > MIN_BLOCK_SIZE)


write.table(q2_synblocks_dist_rfmt[, c('chr_ref', 'start_ref', 'end_ref', 'chr_query', 'start_query', 'end_query', 'orient', 'rID', 'qID')], 
            file = paste0(unlist(strsplit(Q2_BLOCKS, '.tsv')), '_dist_rfmt.txt'),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


draw_linear_laurie(directory = "~/SaFo_paper",
                   paste0('synteny_SaNa_SaFo_SaSaICSASG_', format(Sys.time(), "%Y%m%d.%H%M%S")),
                   sizefile = "~/SaFo_paper/all_chrlengths.txt", 
                   paste0(unlist(strsplit(Q1_BLOCKS, '.tsv')), '_fmt.txt'),
                   paste0(unlist(strsplit(Q2_BLOCKS, '.tsv')), '_dist_rfmt.txt'),
                   #"~/SaFo_paper/querySaSaICSASG_to_refSaFo/maskedSaSaICSASG_on_maskedSaFo_blocks_nucmer_dist_rfmt.txt",
                   colours = rep(c('black', 'blue'), (nrow(all_chrlengths)+1)/2))













SaFo_SaSa_unique_pairs <- unique(Q2_blocks_fmt[, c('chr_ref', 'chr_query')])
SaFo_SaSa_unique_pairs$chr_query <- as.numeric(SaFo_SaSa_unique_pairs$chr_query)
SaFo_SaSa_unique_pairs <- SaFo_SaSa_unique_pairs[order(SaFo_SaSa_unique_pairs$chr_query), ]



write.table(SaFo_SaSa_unique_pairs, file = "~/SaFo_paper/SaFo_SaSa1.3_unique_chr_pairs.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = c('chr_ref_SaFo', 'chr_query_SaSa1.3'))
