# plot synteny with parameters min_dots = 30 and merge_blocks = TRUE

# Test 
Q1_BLOCKS <- "~/SaFo_paper/SaNa_on_SaFo_masked_mindots30_topn1_noover1_merge1_blocks"
Q2_BLOCKS <- "~/SaFo_paper/SaSaICSASG_on_SaFo_masked_mindots30_topn1_noover1_merge1_blocks"

# 1. Import files for reference SaFo --------------------------------------
rchrsizes <- read.delim("~/SaFo_paper/SaFo.chrs_noMt.chrsizes.txt", header=FALSE,
                        col.names = c('chr', 'size'))

rchrs_for_blocks <- read.delim("~/SaFo_paper/SaFo.chrs_for_blocks.txt", header=FALSE,
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



# 2. Import files for query SaNa ------------------------------------------
# Alignment
q1_blocks <- read.delim(Q1_BLOCKS,
                        col.names = c('qChr', 'rChr', 'block', 'qStart', 'qEnd', 'rStart', 'rEnd', 'hits', 'qGenes', 'rGenes', 'qPctGenes', 'rPctGenes'))

for (i in 1:nrow(q1_blocks)){
  q1_blocks$ref_block_size[i] <- (q1_blocks$rEnd[i] - q1_blocks$rStart[i])
  q1_blocks$ref_chr_size[i] <- rchr_lengths$size[rchr_lengths$chrNum == q1_blocks$rChr[i]]
  q1_blocks$ref_chr_prop[i] <- q1_blocks$ref_block_size[i]/q1_blocks$ref_chr_size[i]
  
}

# Add BC from Sutherland et al. 2016
q1_blocks <- merge(x = q1_blocks, y = LG_crrsp, by.x = 'rChr', by.y = 'Claire_assembly', sort = FALSE)

q1_blocks$rChr <- paste0(q1_blocks$rChr, ' (', q1_blocks$Ben_BC, ')')

#q1_blocks <- q1_blocks[order(q1_blocks$rChr),]

# Remove blocks smaller than 1Mb
#q1_blocks <- subset(q1_blocks, ref_block_size > 10000000)

q1_blocks_fmt <- q1_blocks[, c('rChr', 'rStart', 'rEnd', 'qChr', 'qStart', 'qEnd')]
q1_blocks_fmt$orient <- '+'
q1_blocks_fmt$rID <- 'SaFo'
q1_blocks_fmt$qID <- 'SaNa'




write.table(q1_blocks_fmt, file = paste0(Q1_BLOCKS, '_fmt.txt'),
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

q1_chr_order <- c(4,3,2,1,8,5,7,16,21,26,6,9,12,10,15,28,11,20,19,29,31,27,13,23,33,25,24,18,17,35,14,32,34,41,30,36,22,37,38,39,40,42)
q1chr_lengths <- q1chr_lengths[order(match(q1chr_lengths$chrNum, q1_chr_order)),]

#q1chr_lengths <- q1chr_lengths[order(match(q1chr_lengths$chrNum, unique(q1_blocks_fmt$qChr))),]

#all_chrlengths <- rbind(q1chr_lengths, rchr_lengths)
#write.table(all_chrlengths, file = "~/SaFo_paper/querySaNa_to_refSaFo/SaNa_SaFo.chrlengths.txt",
#            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


# 3. Import files for query SaSa ------------------------------------------
# Alignment
q2_blocks <- read.delim(Q2_BLOCKS,
                        col.names = c('qChr', 'rChr', 'block', 'qStart', 'qEnd', 'rStart', 'rEnd', 'hits', 'qGenes', 'rGenes', 'qPctGenes', 'rPctGenes'))

for (i in 1:nrow(q2_blocks)){
  q2_blocks$ref_block_size[i] <- (q2_blocks$rEnd[i] - q2_blocks$rStart[i])
  q2_blocks$ref_chr_size[i] <- rchr_lengths$size[rchr_lengths$chrNum == q2_blocks$rChr[i]]
  q2_blocks$ref_chr_prop[i] <- q2_blocks$ref_block_size[i]/q2_blocks$ref_chr_size[i]
  
}



# Add BC from Sutherland et al. 2016
q2_blocks <- merge(x = q2_blocks, y = LG_crrsp, by.x = 'rChr', by.y = 'Claire_assembly', sort = FALSE)

q2_blocks$rChr <- paste0(q2_blocks$rChr, ' (', q2_blocks$Ben_BC, ')')

#q2_blocks <- subset(q2_blocks, ref_block_size > 10000000)

# Find the biggest block for each chr
#q2_big_blocks <- group_by(q2_blocks, rChr) %>% top_n(1, ref_chr_prop)

#q2_blocks <- rbind(q2_big_blocks, setdiff(q2_blocks, q2_big_blocks))

#q2_blocks <- q2_blocks[order(q2_big_blocks$rChr),]

q2_blocks_fmt <- q2_blocks[, c('rChr', 'rStart', 'rEnd', 'qChr', 'qStart', 'qEnd')]
q2_blocks_fmt$orient <- '+'
q2_blocks_fmt$rID <- 'SaFo'
q2_blocks_fmt$qID <- 'SaSa'


write.table(q2_blocks_fmt, file = paste0(Q2_BLOCKS, '_fmt.txt'),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

#less "synteny/SaFo_SaSa/SaSa.chrs.fasta.fai" | cut -f1,2 | sed -E 's/^(HG.+)\:0\-[0-9]+/\1/'  > "synteny/SaFo_SaSa/SaSa.chrs.chrsizes.txt"
q2chrsizes <- read.delim("C:/Users/Laurie/Documents/SaFo_paper/querySaSaICSASG_to_refSaFo/SaSaICSASG.chrs.chrsizes.txt", header=FALSE,
                         col.names = c('chr', 'size'))

#less "synteny/SaFo_SaSa/SaSa.chrs.corrsp.txt" | sed -E 's/^(HG.+)/:0/-[0-9]+//1/' > "synteny/SaFo_SaSa/SaSa.chrs_for_blocks.txt"
q2chrs_for_blocks <- read.delim("C:/Users/Laurie/Documents/SaFo_paper/querySaSaICSASG_to_refSaFo/SaSaICSASG.chrs_for_blocks.txt", header=FALSE,
                                col.names = c('chrHG', 'chrSymap'))

q2chrs_for_blocks$chrNum <- 1:nrow(q2chrs_for_blocks)


q2chr_lengths <- merge(q2chrsizes, q2chrs_for_blocks, by.x = 'chr', by.y = 'chrSymap', sort = FALSE)
q2chr_lengths <- q2chr_lengths[, c('chrNum', 'size')]
q2chr_lengths$ID <- 'SaSa'

#q2chr_lengths <- q2chr_lengths[order(match(q2chr_lengths$chrNum, unique(q2_big_blocks$qChr))),]
#q2chr_lengths <- q2chr_lengths[order(match(q2chr_lengths$chr, unique(q2_blocks_fmt$qChr))),]

q2_chr_order <- c(1,13,12,24,23,5,16,29,13,10,11,7,16,
                  4,10,1,
                  9, 14, 22, 21,15,1,
                  27,25,20,19,3,6,20,9,28,15,14,9,6,11,18,18,2,17)

q2chr_lengths <- q2chr_lengths[order(match(q2chr_lengths$chrNum, unique(q2_chr_order))),]

# 4. Test syntenyPlotteR --------------------------------------------------
all_chrlengths <- rbind(q1chr_lengths, rchr_lengths, q2chr_lengths)
row.names(all_chrlengths) <- NULL
write.table(all_chrlengths, file = "~/SaFo_paper/all_chrlengths_mindots30.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")




draw_linear_laurie(
  directory = "~/SaFo_paper",
  paste0('symap_SaSaICSASG_SaFo_on_SaFo_masked_mindots30', format(Sys.time(), "%Y%m%d.%H%M%S")),
  sizefile = "~/SaFo_paper/all_chrlengths_mindots30.txt", 
  paste0(Q1_BLOCKS, '_fmt.txt'),
  paste0(Q2_BLOCKS, '_fmt.txt'),
  colours = rep(c('grey20', 'blue'), (nrow(all_chrlengths)+1)/2))



chr_blocks <- (unique(rchr_lengths$chrNum))
cols_blocks <- vector(mode = 'character', length = nrow(rchr_lengths))
#hex_cols <- (viridisLite::viridis(n = length(chr_blocks), option = 'D'))
#hex_cols <- c(viridisLite::viridis(n = length(chr_blocks)/2, option = 'D'), viridisLite::viridis(n = length(chr_blocks)/2, option = 'C'))
col_pair <- rep(c('grey20', 'blue'), (nrow(rchr_lengths)+1)/2)

for (i in 1:length(chr_blocks)) {
  names(cols_blocks)[i] <- chr_blocks[i]
  cols_blocks[i] <- col_pair[i]
}

cols_vec_Q1 <- cols_blocks[match(q1_blocks_fmt$rChr, names(cols_blocks))]
cols_vec_Q2 <- cols_blocks[match(q2_blocks_fmt$rChr, names(cols_blocks))]

full_cols_vec <- c(cols_vec_Q1, cols_blocks, cols_vec_Q2)


draw_linear_laurie(
  directory = "~/SaFo_paper",
  paste0('symap_SaSaICSASG_SaFo_on_SaFo_masked_mindots30', format(Sys.time(), "%Y%m%d.%H%M%S")),
  sizefile = "~/SaFo_paper/all_chrlengths_mindots30.txt", 
  paste0(Q1_BLOCKS, '_fmt.txt'),
  paste0(Q2_BLOCKS, '_fmt.txt'),
  colours = full_cols_vec)



