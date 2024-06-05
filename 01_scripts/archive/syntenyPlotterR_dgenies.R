# 1. Import files for reference SaFo --------------------------------------
rchrsizes <- read.delim("~/SaFo_paper/SaFo.chrs_noMt.chrsizes.txt", header=FALSE,
                        col.names = c('chr', 'size'))

rchrs_for_blocks <- read.delim("~/SaFo_paper/SaFo.chrs_for_blocks.txt", header=FALSE,
                               col.names = c('chrCM', 'chrSymap'))
rchrs_for_blocks$chrNum <- 1:nrow(rchrs_for_blocks)

rchr_lengths <- merge(rchrsizes, rchrs_for_blocks, by.x = 'chr', by.y = 'chrCM', sort = FALSE)
rchr_lengths <- rchr_lengths[, c('chrNum', 'size')]
rchr_lengths$ID <- 'SaFo'

rchr_lengths$chrName <- paste0('Chr', rchr_lengths$chrNum)

# 2. Import files for query SaNa ------------------------------------------
# Alignment
#q1_blocks <- read.delim("~/SaFo_paper/querySaNa_to_refSaFo/results/blocks",
#                        col.names = c('qChr', 'rChr', 'block', 'qStart', 'qEnd', 'rStart', 'rEnd', 'hits', 'qGenes', 'rGenes', 'qPctGenes', 'rPctGenes'))

Q1 <- "~/SaFo_paper/query_SaNa.chrs_noMt_renamed.masked_ref_SaFo.chrs_noMt_renamed.masked_assoc.tsv" 
q1_blocks <- read.delim(Q1)


#q1_blocks$Q.Chr <- paste0('Chr', q1_blocks$Query)
#q1_blocks$R.Chr <- paste0('Chr', q1_blocks$Target)

for (i in 1:nrow(q1_blocks)){
  q1_blocks$ref_block_size[i] <- (q1_blocks$T.stop[i] - q1_blocks$T.start[i])
  q1_blocks$ref_chr_size[i] <- rchr_lengths$size[rchr_lengths$chrName == q1_blocks$Target[i]]
  q1_blocks$ref_chr_prop[i] <- q1_blocks$ref_block_size[i]/q1_blocks$ref_chr_size[i]
  
}

q1_blocks <- q1_blocks[order(q1_blocks$Target),]

# Remove blocks smaller than 100 kb
q1_blocks <- subset(q1_blocks, ref_block_size > 5000000)

#q1_blocks_fmt <- q1_blocks[, c('rChr', 'rStart', 'rEnd', 'qChr', 'qStart', 'qEnd')]
q1_blocks_fmt <- q1_blocks[, c('Target', 'T.start', 'T.stop', 'Query', 'Q.start', 'Q.stop')]
q1_blocks_fmt$orient <- '+'
q1_blocks_fmt$rID <- 'SaFo'
q1_blocks_fmt$qID <- 'SaNa'


# Remove leading Chr
q1_blocks_fmt$Target <- sapply(X = q1_blocks_fmt$Target, FUN = function(x) substr(x, start = 4, stop = 6))
q1_blocks_fmt$Query <- sapply(X = q1_blocks_fmt$Query, FUN = function(x) substr(x, start = 4, stop = 6))


write.table(q1_blocks_fmt, file = "~/SaFo_paper/dgenies_SaNa_on_SaFo.masked.blocks.txt",
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


#q1chr_lengths <- q1chr_lengths[order(match(q1chr_lengths$chrNum, unique(q1_blocks_fmt$Query))),]
q1_chr_order <- c(4,3,2,1,8,5,7,16,21,26,6,9,12,10,15,28,11,20,19,29,31,27,13,23,33,25,24,18,17,35,14,32,34,41,30,36,22,37,38,39,40,42)
q1chr_lengths <- q1chr_lengths[order(match(q1chr_lengths$chrNum, q1_chr_order)),]


# 2. Import files for query SaSa ------------------------------------------
# Alignment
#q1_blocks <- read.delim("~/SaFo_paper/querySaNa_to_refSaFo/results/blocks",
#                        col.names = c('qChr', 'rChr', 'block', 'qStart', 'qEnd', 'rStart', 'rEnd', 'hits', 'qGenes', 'rGenes', 'qPctGenes', 'rPctGenes'))

Q2 <- "~/SaFo_paper/query_SaSa.chrs_renamed.masked_ref_SaFo.chrs_noMt_renamed.masked_assoc.tsv" 
q2_blocks <- read.delim(Q2)

for (i in 1:nrow(q2_blocks)){
  q2_blocks$ref_block_size[i] <- (q2_blocks$T.stop[i] - q2_blocks$T.start[i])
  #q2_blocks$ref_chr_size[i] <- rchr_lengths$size[rchr_lengths$chrName == q2_blocks$Target[i]]
  q2_blocks$ref_chr_size[i] <- rchr_lengths$size[rchr_lengths$chrName == q2_blocks$Target[i]]
  q2_blocks$ref_chr_prop[i] <- q2_blocks$ref_block_size[i]/q2_blocks$ref_chr_size[i]
  
}

q2_blocks <- q2_blocks[order(q2_blocks$Target),]

# Remove blocks smaller than 100 kb
q2_blocks <- subset(q2_blocks, ref_block_size > 5000000)

#q2_blocks_fmt <- q2_blocks[, c('rChr', 'rStart', 'rEnd', 'qChr', 'qStart', 'qEnd')]
q2_blocks_fmt <- q2_blocks[, c('Target', 'T.start', 'T.stop', 'Query', 'Q.start', 'Q.stop')]
q2_blocks_fmt$orient <- '+'
q2_blocks_fmt$rID <- 'SaFo'
q2_blocks_fmt$qID <- 'SaSa'

# Remove leading Chr
q2_blocks_fmt$Target <- sapply(X = q2_blocks_fmt$Target, FUN = function(x) substr(x, start = 4, stop = 6))
q2_blocks_fmt$Query <- sapply(X = q2_blocks_fmt$Query, FUN = function(x) substr(x, start = 4, stop = 6))



write.table(q2_blocks_fmt, file = "~/SaFo_paper/dgenies_SaSa_on_SaFo.masked.blocks.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

#less "synteny/SaFo_SaSa/SaSa.chrs_noMt.fasta.fai" | cut -f1,2 | sed -E 's/^(CM.+)\:0\-[0-9]+/\1/'  > "synteny/SaFo_SaSa/SaSa.chrs_noMt.chrsizes.txt"
q2chrsizes <- read.delim("~/SaFo_paper/querySaSa_to_refSaFo/SaSa.chrs.chrsizes.txt", header=FALSE,
                         col.names = c('chr', 'size'))

#less "synteny/SaFo_SaSa/SaSa.chrs_noMt.corrsp.txt" | sed -E 's/^(CM.+)\:0\-[0-9]+/\1/' > "synteny/SaFo_SaSa/SaSa.chrs_for_blocks.txt"
q2chrs_for_blocks <- read.delim("~/SaFo_paper/querySaSa_to_refSaFo/SaSa.chrs_for_blocks.txt", header=FALSE,
                                col.names = c('chrCM', 'chrSymap'))

q2chrs_for_blocks$chrNum <- 1:nrow(q2chrs_for_blocks)


q2chr_lengths <- merge(q2chrsizes, q2chrs_for_blocks, by.x = 'chr', by.y = 'chrCM', sort = FALSE)
q2chr_lengths <- q2chr_lengths[, c('chrNum', 'size')]
q2chr_lengths$ID <- 'SaSa'


#q2chr_lengths <- q2chr_lengths[order(match(q2chr_lengths$chrNum, unique(q2_blocks_fmt$Query))),]

#q2_chr_order <- c(1,13,12,24,23,5,16,29,13,10,11,7,16,
#                  4,10,1,
#                  9, 14, 22, 21,15,1,
#                  27,25,20,19,3,6,20,9,28,15,14,9,6,11,18,18,2,17)

#q2chr_lengths <- q2chr_lengths[order(match(q2chr_lengths$chrNum, unique(q2_chr_order))),]
q2chr_lengths <- q2chr_lengths[order(match(q2chr_lengths$chrNum, 
                                           q2_blocks_fmt$Query[order(as.numeric(q2_blocks_fmt$Target))]
                                           )), ]


# 4. Test syntenyPlotteR --------------------------------------------------

all_chrlengths <- rbind(q1chr_lengths, rchr_lengths[, c('chrNum', 'size', 'ID')], q2chr_lengths)

write.table(all_chrlengths, file = "~/SaFo_paper/all_chrlengths_dgenies.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")



library(syntenyPlotteR)
draw.linear(directory = "~/SaFo_paper",
            paste0('test', format(Sys.time(), "%Y%m%d.%H%M%S")),
            sizefile = "~/SaFo_paper/all_chrlengths_dgenies.txt", 
            "~/SaFo_paper/dgenies_SaNa_on_SaFo.masked.blocks.txt",
            "~/SaFo_paper/dgenies_SaSa_on_SaFo.masked.blocks.txt",
            colours = rep(c('black', 'blue'), (nrow(all_chrlengths)+1)/2))

draw_linear_laurie(
  directory = "~/SaFo_paper",
  paste0('test', format(Sys.time(), "%Y%m%d.%H%M%S")),
  sizefile = "~/SaFo_paper/all_chrlengths.txt", 
  "~/SaFo_paper/querySaNa_to_refSaFo/results/blocks_fmt.txt",
  "~/SaFo_paper/querySaSa_to_refSaFo/results/blocks_fmt.txt",
  colours = rep(c('grey20', 'blue'), (nrow(all_chrlengths)+1)/2))

