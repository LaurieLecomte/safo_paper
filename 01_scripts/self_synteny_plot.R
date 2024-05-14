# Plot self synteny from Symap output

CHR_SIZES <- "~/SaFo_paper/SaFo.chrs_noMt.chrsizes.txt"
CHR_BLOCKS <- "~/SaFo_paper/SaFo.chrs_for_blocks.txt"

SELF_BLOCKS <- "~/SaFo_paper/maskedSaFo_maskedSaFo_blocks"

MIN_BLOCK_SIZE <- 1000000

# 1. Import files for reference SaFo --------------------------------------
rchrsizes <- read.delim(CHR_SIZES, header=FALSE,
                        col.names = c('chr', 'size'))

rchrs_for_blocks <- read.delim(CHR_BLOCKS, header=FALSE,
                               col.names = c('chrCM', 'chrSymap'))
rchrs_for_blocks$chrNum <- 1:nrow(rchrs_for_blocks)

rchr_lengths <- merge(rchrsizes, rchrs_for_blocks, by.x = 'chr', by.y = 'chrCM', sort = FALSE)
rchr_lengths <- rchr_lengths[, c('chrNum', 'size')]
rchr_lengths$ID <- 'SaFo'


# 2. Format input blocks from symap -----------------------------------

self_syn_blocks <- read.delim(SELF_BLOCKS)
self_syn_blocks$block_size1 <- self_syn_blocks$end1 - self_syn_blocks$start1
self_syn_blocks$block_size2 <- self_syn_blocks$end2 - self_syn_blocks$start2


self_syn_blocks_filt <- subset(self_syn_blocks, block_size1 > MIN_BLOCK_SIZE & block_size2 > MIN_BLOCK_SIZE)

bed1 <- self_syn_blocks_filt[, c('grp1', 'start1', 'end1')]

bed2 <- self_syn_blocks_filt[, c('grp2', 'start2', 'end2')]




# 3. Plot self synteny using circos ---------------------------------------
library(circlize)
# Initialize plot space and chr bins
circos.clear()
col_text <- "grey40"
circos.par("track.height"=0.8, gap.degree=1, cell.padding=c(0, 0, 0, 0))

rchr_lengths$start <- 0
rchr_lengths <- rchr_lengths[, c('chrNum', 'start', 'size')]

circos.initialize(factors=rchr_lengths$chrNum, 
                  xlim= rchr_lengths[, c('start', 'size')])

circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr, cex=0.5, col=col_text, 
              facing="bending.inside", niceFacing=TRUE)
}, bg.col="grey90", bg.border=F, track.height=0.06)


# Create color scheme, 1 for each grp1 chr
chr_blocks <- sort(unique(self_syn_blocks_filt$grp1)) 
# Assign a color to each svtype in a named vector
cols_blocks <- vector(mode = 'character', length = length(chr_blocks))
#hex_cols <- (viridisLite::viridis(n = length(chr_blocks), option = 'D'))
hex_cols <- c(viridisLite::viridis(n = length(chr_blocks)/2, option = 'D'), viridisLite::viridis(n = length(chr_blocks)/2, option = 'C'))
for (i in 1:length(chr_blocks)) {
  names(cols_blocks)[i] <- chr_blocks[i]
  cols_blocks[i] <- hex_cols[i]
}

cols_vec <- cols_blocks[match(self_syn_blocks_filt$grp1, names(cols_blocks))]
cols_vec <- scales::alpha(cols_vec, alpha = 0.75)

circos.genomicLink(bed1, bed2, col = cols_vec)




# show biggest blocks first, then add smaller ones on top
self_syn_blocks_filt_large <- subset(self_syn_blocks_filt, block_size1 > 5000000 & block_size2 > 5000000)
bed1_large <- self_syn_blocks_filt_large[, c('grp1', 'start1', 'end1')]
bed2_large <- self_syn_blocks_filt_large[, c('grp2', 'start2', 'end2')]


# Create color scheme, 1 for each grp1 chr
chr_blocks_large <- sort(unique(self_syn_blocks_filt_large$grp1)) 
# Assign a color to each svtype in a named vector
cols_blocks_large <- vector(mode = 'character', length = length(chr_blocks_large))
hex_cols_large <- (viridisLite::viridis(n = length(chr_blocks_large), option = 'D'))
for (i in 1:length(chr_blocks_large)) {
  names(cols_blocks_large)[i] <- chr_blocks_large[i]
  cols_blocks_large[i] <- hex_cols_large[i]
}

cols_vec_large <- cols_blocks_large[match(bed1_large$grp1, names(cols_blocks_large))]
cols_vec_large <- scales::alpha(cols_vec_large, alpha = 0.75)



self_syn_blocks_filt_small <- subset(self_syn_blocks_filt, block_size1 <= 5000000 & block_size2 <= 5000000)
bed1_small <- self_syn_blocks_filt_small[, c('grp1', 'start1', 'end1')]
bed2_small <- self_syn_blocks_filt_small[, c('grp2', 'start2', 'end2')]

# Create color scheme, 1 for each grp1 chr
chr_blocks_small <- sort(unique(self_syn_blocks_filt_small$grp1)) 
# Assign a color to each svtype in a named vector
cols_blocks_small <- vector(mode = 'character', length = length(chr_blocks_small))
hex_cols_small <- (viridisLite::viridis(n = length(chr_blocks_small), option = 'C'))
for (i in 1:length(chr_blocks_small)) {
  names(cols_blocks_small)[i] <- chr_blocks_small[i]
  cols_blocks_small[i] <- hex_cols_small[i]
}

cols_vec_small <- cols_blocks_small[match(bed1_small$grp1, names(cols_blocks_small))]
cols_vec_small <- scales::alpha(cols_vec_small, alpha = 0.75)





circos.clear()
col_text <- "grey40"
circos.par("track.height"=0.8, gap.degree=1, cell.padding=c(0, 0, 0, 0))

rchr_lengths$start <- 0
rchr_lengths <- rchr_lengths[, c('chrNum', 'start', 'size')]

circos.initialize(factors=rchr_lengths$chrNum, 
                  xlim= rchr_lengths[, c('start', 'size')])

circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr, cex=0.5, col=col_text, 
              facing="bending.inside", niceFacing=TRUE)
}, bg.col="grey90", bg.border=F, track.height=0.06)


circos.genomicLink(bed1_large, bed2_large, 
                   col = cols_vec_large, lwd = 0.05, border = 'grey70')
circos.genomicLink(bed1_small, bed2_small, 
                   col = cols_vec_small)












# Create color scheme, 1 for each grp1 chr
chr_blocks <- sort(unique(self_syn_blocks_filt_large$grp1)) 
# Assign a color to each svtype in a named vector
cols_blocks <- vector(mode = 'character', length = length(chr_blocks))
hex_cols <- (viridisLite::viridis(n = length(chr_blocks), option = 'D'))
for (i in 1:length(chr_blocks)) {
  names(cols_blocks)[i] <- chr_blocks[i]
  cols_blocks[i] <- hex_cols[i]
}

cols_vec <- cols_blocks[match(bed1_large$grp1, names(cols_blocks))]
cols_vec <- scales::alpha(cols_vec, alpha = 0.5)




