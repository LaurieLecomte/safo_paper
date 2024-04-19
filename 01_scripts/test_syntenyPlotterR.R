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
q1_blocks <- read.delim("~/SaFo_paper/querySaNa_to_refSaFo/results/blocks",
                     col.names = c('qChr', 'rChr', 'block', 'qStart', 'qEnd', 'rStart', 'rEnd', 'hits', 'qGenes', 'rGenes', 'qPctGenes', 'rPctGenes'))


for (i in 1:nrow(q1_blocks)){
  q1_blocks$ref_block_size[i] <- (q1_blocks$rEnd[i] - q1_blocks$rStart[i])
  q1_blocks$ref_chr_size[i] <- rchr_lengths$size[rchr_lengths$chrNum == q1_blocks$rChr[i]]
  q1_blocks$ref_chr_prop[i] <- q1_blocks$ref_block_size[i]/q1_blocks$ref_chr_size[i]
  
}

q1_blocks <- q1_blocks[order(q1_blocks$rChr),]

# Remove blocks smaller than 100 kb
q1_blocks <- subset(q1_blocks, ref_block_size > 5000000)

q1_blocks_fmt <- q1_blocks[, c('rChr', 'rStart', 'rEnd', 'qChr', 'qStart', 'qEnd')]
q1_blocks_fmt$orient <- '+'
q1_blocks_fmt$rID <- 'SaFo'
q1_blocks_fmt$qID <- 'SaNa'



write.table(q1_blocks_fmt, file = "~/SaFo_paper/querySaNa_to_refSaFo/results/blocks_fmt.txt",
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


q1chr_lengths <- q1chr_lengths[order(match(q1chr_lengths$chrNum, unique(q1_blocks_fmt$qChr))),]

#all_chrlengths <- rbind(q1chr_lengths, rchr_lengths)
#write.table(all_chrlengths, file = "~/SaFo_paper/querySaNa_to_refSaFo/SaNa_SaFo.chrlengths.txt",
#            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


# 3. Import files for query SaSa ------------------------------------------
# Alignment
q2_blocks <- read.delim("~/SaFo_paper/querySaSa_to_refSaFo/results/blocks",
                     col.names = c('qChr', 'rChr', 'block', 'qStart', 'qEnd', 'rStart', 'rEnd', 'hits', 'qGenes', 'rGenes', 'qPctGenes', 'rPctGenes'))

for (i in 1:nrow(q2_blocks)){
  q2_blocks$ref_block_size[i] <- (q2_blocks$rEnd[i] - q2_blocks$rStart[i])
  q2_blocks$ref_chr_size[i] <- rchr_lengths$size[rchr_lengths$chrNum == q2_blocks$rChr[i]]
  q2_blocks$ref_chr_prop[i] <- q2_blocks$ref_block_size[i]/q2_blocks$ref_chr_size[i]
  
}




q2_blocks <- subset(q2_blocks, ref_block_size > 5000000)

# Find the biggest block for each chr
q2_big_blocks <- group_by(q2_blocks, rChr) %>% top_n(1, ref_chr_prop)

q2_blocks <- rbind(q2_big_blocks, setdiff(q2_blocks, q2_big_blocks))

q2_blocks <- q2_blocks[order(q2_big_blocks$rChr),]

q2_blocks_fmt <- q2_blocks[, c('rChr', 'rStart', 'rEnd', 'qChr', 'qStart', 'qEnd')]
q2_blocks_fmt$orient <- '+'
q2_blocks_fmt$rID <- 'SaFo'
q2_blocks_fmt$qID <- 'SaSa'


write.table(q2_blocks_fmt, file = "~/SaFo_paper/querySaSa_to_refSaFo/results/blocks_fmt.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

#less "synteny/SaFo_SaSa/SaSa.chrs.fasta.fai" | cut -f1,2 | sed -E 's/^(HG.+)\:0\-[0-9]+/\1/'  > "synteny/SaFo_SaSa/SaSa.chrs.chrsizes.txt"
q2chrsizes <- read.delim("~/SaFo_paper/querySaSa_to_refSaFo/SaSa.chrs.chrsizes.txt", header=FALSE,
                         col.names = c('chr', 'size'))

#less "synteny/SaFo_SaSa/SaSa.chrs.corrsp.txt" | sed -E 's/^(HG.+)\:0\-[0-9]+/\1/' > "synteny/SaFo_SaSa/SaSa.chrs_for_blocks.txt"
q2chrs_for_blocks <- read.delim("~/SaFo_paper/querySaSa_to_refSaFo/SaSa.chrs_for_blocks.txt", header=FALSE,
                                col.names = c('chrHG', 'chrSymap'))

q2chrs_for_blocks$chrNum <- 1:nrow(q2chrs_for_blocks)


q2chr_lengths <- merge(q2chrsizes, q2chrs_for_blocks, by.x = 'chr', by.y = 'chrHG', sort = FALSE)
q2chr_lengths <- q2chr_lengths[, c('chrNum', 'size')]
q2chr_lengths$ID <- 'SaSa'

#q2chr_lengths <- q2chr_lengths[order(match(q2chr_lengths$chrNum, unique(q2_big_blocks$qChr))),]
q2chr_lengths <- q2chr_lengths[order(match(q2chr_lengths$chr, unique(q2_blocks_fmt$qChr))),]

# 4. Test syntenyPlotteR --------------------------------------------------
all_chrlengths <- rbind(q1chr_lengths, rchr_lengths, q2chr_lengths)

write.table(all_chrlengths, file = "~/SaFo_paper/all_chrlengths.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")



library(syntenyPlotteR)
draw.linear(directory = "~/SaFo_paper",
            paste0('test', format(Sys.time(), "%Y%m%d.%H%M%S")),
            sizefile = "~/SaFo_paper/all_chrlengths.txt", 
            "~/SaFo_paper/querySaNa_to_refSaFo/results/blocks_fmt.txt",
            "~/SaFo_paper/querySaSa_to_refSaFo/results/blocks_fmt.txt",
            colours = rep(c('black', 'blue'), (nrow(all_chrlengths)+1)/2))

draw_linear_laurie(
  directory = "~/SaFo_paper",
  paste0('test', format(Sys.time(), "%Y%m%d.%H%M%S")),
            sizefile = "~/SaFo_paper/all_chrlengths.txt", 
            "~/SaFo_paper/querySaNa_to_refSaFo/results/blocks_fmt.txt",
            "~/SaFo_paper/querySaSa_to_refSaFo/results/blocks_fmt.txt",
            colours = rep(c('grey20', 'blue'), (nrow(all_chrlengths)+1)/2))












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

draw_linear_laurie(directory = "~/SaFo_paper",
            'test',
            sizefile = "~/SaFo_paper/all_chrlengths.txt", 
            "~/SaFo_paper/querySaNa_to_refSaFo/results/blocks_fmt.txt",
            "~/SaFo_paper/querySaMa_to_refSaFo/results/blocks_fmt.txt",
            colours = rep(c('black', 'blue'), (nrow(all_chrlengths)+1)/2))














# All sizes together ------------------------------------------------------
rchrsizes <- read.delim("~/SaFo_paper/SaFo.chrs_noMt.chrsizes.txt", header=FALSE,
                        col.names = c('chr', 'size'))

rchrs_for_blocks <- read.delim("~/SaFo_paper/SaFo.chrs_for_blocks.txt", header=FALSE,
                               col.names = c('chrCM', 'chrSymap'))
rchrs_for_blocks$chrNum <- 1:nrow(rchrs_for_blocks)

rchr_lengths <- merge(rchrsizes, rchrs_for_blocks, by.x = 'chr', by.y = 'chrCM', sort = FALSE)
rchr_lengths$ID <- 'SaFo'

#less "synteny/SaFo_SaNa/SaNa.chrs_noMt.fasta.fai" | cut -f1,2 | sed -E 's/^(CM.+)\:0\-[0-9]+/\1/'  > "synteny/SaFo_SaNa/SaNa.chrs_noMt.chrsizes.txt"
q1chrsizes <- read.delim("~/SaFo_paper/querySaNa_to_refSaFo/SaNa.chrs_noMt.chrsizes.txt", header=FALSE,
                         col.names = c('chr', 'size'))

#less "synteny/SaFo_SaNa/SaNa.chrs_noMt.corrsp.txt" | sed -E 's/^(CM.+)\:0\-[0-9]+/\1/' > "synteny/SaFo_SaNa/SaNa.chrs_for_blocks.txt"
q1chrs_for_blocks <- read.delim("~/SaFo_paper/querySaNa_to_refSaFo/SaNa.chrs_for_blocks.txt", header=FALSE,
                                col.names = c('chr', 'chrSymap'))

q1chrs_for_blocks$chrNum <- 1:nrow(q1chrs_for_blocks)


q1chr_lengths <- merge(q1chrsizes, q1chrs_for_blocks, by = 'chr', sort = FALSE)
#q1chr_lengths <- q1chr_lengths[, c('chrNum', 'size')]
q1chr_lengths$ID <- 'SaNa'


#q1chr_lengths <- q1chr_lengths[order(match(q1chr_lengths$chrNum, unique(q1_blocks_fmt$qChr))),]



#less "synteny/SaFo_SaSa/SaSa.chrs.fasta.fai" | cut -f1,2 | sed -E 's/^(HG.+)\:0\-[0-9]+/\1/'  > "synteny/SaFo_SaSa/SaSa.chrs.chrsizes.txt"
q2chrsizes <- read.delim("~/SaFo_paper/querySaSa_to_refSaFo/SaSa.chrs.chrsizes.txt", header=FALSE,
                         col.names = c('chr', 'size'))

#less "synteny/SaFo_SaSa/SaSa.chrs.corrsp.txt" | sed -E 's/^(HG.+)\:0\-[0-9]+/\1/' > "synteny/SaFo_SaSa/SaSa.chrs_for_blocks.txt"
q2chrs_for_blocks <- read.delim("~/SaFo_paper/querySaSa_to_refSaFo/SaSa.chrs_for_blocks.txt", header=FALSE,
                                col.names = c('chr', 'chrSymap'))

q2chrs_for_blocks$chrNum <- 1:nrow(q2chrs_for_blocks)


q2chr_lengths <- merge(q2chrsizes, q2chrs_for_blocks, by = 'chr', sort = FALSE)
#q2chr_lengths <- q2chr_lengths[, c('chrNum', 'size')]
q2chr_lengths$ID <- 'SaSa'

#q2chr_lengths <- q2chr_lengths[order(match(q2chr_lengths$chrNum, unique(q2_big_blocks$qChr))),]


#all_lengths_table <- rbind(q1chr_lengths, rchr_lengths, q2chr_lengths)



# Query 1
q1_blocks <- read.delim("~/SaFo_paper/querySaNa_to_refSaFo/results/blocks",
                        col.names = c('qChr', 'rChr', 'block', 'qStart', 'qEnd', 'rStart', 'rEnd', 'hits', 'qGenes', 'rGenes', 'qPctGenes', 'rPctGenes'))

q1_blocks$qID <- 'SaNa'


for (i in 1:nrow(q1_blocks)){
  q1_blocks$ref_block_size[i] <- (q1_blocks$rEnd[i] - q1_blocks$rStart[i])
  q1_blocks$ref_chr_size[i] <- rchr_lengths$size[rchr_lengths$chrNum == q1_blocks$rChr[i]]
  q1_blocks$ref_chr_prop[i] <- q1_blocks$ref_block_size[i]/q1_blocks$ref_chr_size[i]
  q1_blocks$rChrCM[i] <- rchr_lengths$chr[rchr_lengths$chrNum == q1_blocks$rChr[i]]
  q1_blocks$qChrCM[i] <- q1chr_lengths$chr[q1chr_lengths$chrNum == q1_blocks$qChr[i]]
  
}

q1_blocks <- q1_blocks[order(q1_blocks$rChr),]

# Remove blocks smaller than 100 kb
q1_blocks <- subset(q1_blocks, ref_block_size > 20000000)

#q1_blocks_fmt <- q1_blocks[, c('rChrCM', 'rStart', 'rEnd', 'qChrCM', 'qStart', 'qEnd')]
q1_blocks_fmt <- q1_blocks[, c('rChr', 'rStart', 'rEnd', 'qChr', 'qStart', 'qEnd')]

q1_blocks_fmt$orient <- '+'
q1_blocks_fmt$rID <- 'SaFo'
q1_blocks_fmt$qID <- 'SaNa'


write.table(q1_blocks_fmt, file = "~/SaFo_paper/querySaNa_to_refSaFo/results/blocks_fmt1.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


#q1chr_lengths <- q1chr_lengths[order(match(q1chr_lengths$chr, unique(q1_blocks_fmt$qChrCM))),]
q1chr_lengths <- q1chr_lengths[order(match(q1chr_lengths$chr, unique(q1_blocks_fmt$qChr))),]





# Query 2
q2_blocks <- read.delim("~/SaFo_paper/querySaSa_to_refSaFo/results/blocks",
                        col.names = c('qChr', 'rChr', 'block', 'qStart', 'qEnd', 'rStart', 'rEnd', 'hits', 'qGenes', 'rGenes', 'qPctGenes', 'rPctGenes'))

q2_blocks$qID <- 'SaSa'



for (i in 1:nrow(q2_blocks)){
  q2_blocks$ref_block_size[i] <- (q2_blocks$rEnd[i] - q2_blocks$rStart[i])
  q2_blocks$ref_chr_size[i] <- rchr_lengths$size[rchr_lengths$chrNum == q2_blocks$rChr[i]]
  q2_blocks$ref_chr_prop[i] <- q2_blocks$ref_block_size[i]/q2_blocks$ref_chr_size[i]
  q2_blocks$rChrCM[i] <- rchr_lengths$chr[rchr_lengths$chrNum == q2_blocks$rChr[i]]
  q2_blocks$qChrCM[i] <- q2chr_lengths$chr[q2chr_lengths$chrNum == q2_blocks$qChr[i]]
  
}

q2_blocks <- q2_blocks[order(q2_blocks$rChr),]

# Remove blocks smaller than 100 kb
q2_blocks <- subset(q2_blocks, ref_block_size > 20000000)


# Find the biggest block for each chr
q2_big_blocks <- group_by(q2_blocks, rChr) %>% top_n(1, ref_chr_prop)

#q2_big_blocks <- group_by(q2_blocks, rChr) %>% summarise(main_chr = max(ref_chr_pro), qChr = qChr)

q2_blocks <- rbind(q2_big_blocks, setdiff(q2_blocks, q2_big_blocks))

#q2_blocks <- q2_blocks[order(q2_big_blocks$rChr),]

#q2_blocks_fmt <- q2_blocks[, c('rChrCM', 'rStart', 'rEnd', 'qChrCM', 'qStart', 'qEnd')]
q2_blocks_fmt <- q2_blocks[, c('rChr', 'rStart', 'rEnd', 'qChr', 'qStart', 'qEnd')]

q2_blocks_fmt$orient <- '+'
q2_blocks_fmt$rID <- 'SaFo'
q2_blocks_fmt$qID <- 'SaSa'



write.table(q2_blocks_fmt, file = "~/SaFo_paper/querySaSa_to_refSaFo/results/blocks_fmt1.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


#q2chr_lengths <- q2chr_lengths[order(match(q2chr_lengths$chr, unique(q2_blocks_fmt$qChrCM))),]
q2chr_lengths <- q2chr_lengths[order(match(q2chr_lengths$chr, unique(q2_blocks_fmt$qChr))),]





all_chrlengths <- rbind(q1chr_lengths, 
                        rchr_lengths, 
                        q2chr_lengths)

write.table(all_chrlengths[, c('chr', 'size', 'ID')], file = "~/SaFo_paper/all_chrlengths_CM.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")



library(syntenyPlotteR)
draw.linear(directory = "~/SaFo_paper",
            paste0('test', format(Sys.time(), "%Y%m%d.%H%M%S")),
            sizefile = "~/SaFo_paper/all_chrlengths_CM.txt", 
            "~/SaFo_paper/querySaNa_to_refSaFo/results/blocks_fmt1.txt",
            "~/SaFo_paper/querySaSa_to_refSaFo/results/blocks_fmt1.txt",
            colours = rep(c('black', 'blue'), (nrow(all_chrlengths)+1)/2))



draw_linear_laurie(directory = "~/SaFo_paper",
            paste0('test', format(Sys.time(), "%Y%m%d.%H%M%S")),
            sizefile = "~/SaFo_paper/all_chrlengths_CM.txt", 
            "~/SaFo_paper/querySaNa_to_refSaFo/results/blocks_fmt1.txt",
            "~/SaFo_paper/querySaSa_to_refSaFo/results/blocks_fmt1.txt",
            colours = rep(c('black', 'blue'), (nrow(all_chrlengths)+1)/2))






# With manual order -------------------------------------------------------
q1_chr_order <- c(4,3,2,1,8,5,7,16,21,26,6,9,12,10,15,28,11,20,19,29,31,27,13,23,33,25,24,18,17,35,14,32,34,41,30,36,22,37,38,39,40,42)




