# Plot self synteny from Symap output

library(tidyr)
library(dplyr)
library(circlize)

options(scipen = 999)

# Chromosomes
CHR_SIZES <- "~/SaFo_paper/SaFo.chrs_noMt.chrsizes.txt"
CHR_BLOCKS <- "~/SaFo_paper/SaFo.chrs_for_blocks.txt"

# Synteny blocks
#SELF_BLOCKS <- "~/SaFo_paper/maskedSaFo_maskedSaFo_blocks"
SELF_BLOCKS <- "~/SaFo_paper/SaFo_self_masked_mindots30_topn2_blocks"
#SELF_BLOCKS <- "~/SaFo_paper/maskedSaFo_on_maskedSaFo_mindots30_blocks"

MIN_BLOCK_SIZE <- 1000000 # blocks smaller than this will not be shown in the plot

# Depth/coverage
#DEPTH_BED <- "~/SaFo_paper/modepth_1Mb.txt.regions.bed.gz"
DEPTH_BED <- "~/SaFo_paper/mosdepth_SaFo_cat_1_1Mb.regions.bed.gz"

# Identity/similarity
IDY_BED <- "~/SaFo_paper/homolog_blocks_identity_0_win1Mb.bed"


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

## Add block size
self_syn_blocks$block_size1 <- self_syn_blocks$end1 - self_syn_blocks$start1
self_syn_blocks$block_size2 <- self_syn_blocks$end2 - self_syn_blocks$start2

## Get prop of chromosome base pairs located in a homeolog region
TOT_CHR_SIZE <- 2218585108
sum(self_syn_blocks$block_size1) / TOT_CHR_SIZE

# Remove blocks smaller than min block size
self_syn_blocks_filt <- subset(self_syn_blocks, block_size1 > MIN_BLOCK_SIZE & block_size2 > MIN_BLOCK_SIZE)

# Split into large and small blocks
self_syn_blocks_filt_large <- subset(self_syn_blocks_filt, block_size1 > 5000000 & block_size2 > 5000000)
self_syn_blocks_filt_small <- subset(self_syn_blocks_filt, block_size1 <= 5000000 & block_size2 <= 5000000)

bed1_large <- self_syn_blocks_filt_large[, c('grp1', 'start1', 'end1')]
bed2_large <- self_syn_blocks_filt_large[, c('grp2', 'start2', 'end2')]

bed1_small <- self_syn_blocks_filt_small[, c('grp1', 'start1', 'end1')]
bed2_small <- self_syn_blocks_filt_small[, c('grp2', 'start2', 'end2')]



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





# 3. Import and format coverage/depth info ---------------------------
depth <- read.delim(DEPTH_BED, header=FALSE,
                    col.names = c('CHROM', 'START', 'STOP', 'MEAN_DEPTH'))
depth$CHROM <- substr(depth$CHROM, start = 4, stop = 6)

blue_no <- adjustcolor( "blue", alpha.f = 0.3)
red_yes <- adjustcolor( "red", alpha.f = 0.3)


depth$too_high <- ifelse(depth$MEAN_DEPTH > mean(depth$MEAN_DEPTH) + 2*(sd(depth$MEAN_DEPTH)),
                         yes = red_yes , no = blue_no)

# 4. Import and format similarity info ------------------------------------
identity_win <- read.delim(IDY_BED, header=FALSE,
                           col.names = c('CHROM', 'START', 'STOP', 'IDY', 'IDY_W'))
identity_win$MID <- identity_win$START + ((identity_win$STOP - identity_win$START)/2)
identity_win$BP <- identity_win$STOP -identity_win$START




# 5. Plot everything together ---------------------------------------------

# Initialize circos plot
circos.clear()
circos.par("track.height"=0.8, gap.degree=0.5, cell.padding=c(0, 0, 0, 0))

rchr_lengths$start <- 0
rchr_lengths <- rchr_lengths[, c('chrNum', 'start', 'size')]

circos.initialize(factors=rchr_lengths$chrNum, 
                  xlim= rchr_lengths[, c('start', 'size')])

# Add chromosome sectors
circos.track(ylim=c(0, 1), 
             panel.fun=function(x, y) {
               chr=CELL_META$sector.index
               xlim=CELL_META$xlim
               ylim=CELL_META$ylim
               circos.text(mean(xlim), mean(ylim), chr, cex=0.5, col='grey10', 
                           facing="bending.inside", niceFacing=TRUE)
             }, bg.col="grey60", bg.border=F, track.height=0.04)

# Add track for sequencing depth
circos.genomicTrackPlotRegion(depth, 
                              ylim = c(min(depth$MEAN_DEPTH), max(depth$MEAN_DEPTH + 2)),
                              panel.fun = function(region, value, ...) {
                                #col = value$too_high
                                col = value$too_high
                                cex = sqrt(value$MEAN_DEPTH)/30
                                circos.genomicPoints(region, value, 
                                                   #ybottom = 0, 
                                                   #ytop = 500, 
                                                   col = col, 
                                                   cex = cex, pch = 16,
                                                   border = NA)
                           
                              }, bg.border= 'black',
                              track.height=0.06)

# Add track for %similarity
##Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('white', 'grey80',  'grey50', 'black'))

identity_win$cols_idy <- rbPal(20)[as.numeric(cut(identity_win$IDY_W, breaks = 20))]

## Add to plot
circos.genomicTrackPlotRegion(identity_win, ylim = c(0, 1),
                              panel.fun = function(region, value, ...) {
                                col = value$cols_idy
                                circos.genomicRect(region, value, 
                                                   ybottom = 0, ytop = 1, 
                                                   col = col, border = NA)
                              }, bg.border= 'black', 
                              bg.col = 'steelblue1',
                              track.height=0.06)


# Add inner links for homeologs
circos.genomicLink(bed1_large, bed2_large, 
                   col = cols_vec_large, 
                   #lwd = 0.04, 
                   #border = 'grey70'
                   )
circos.genomicLink(bed1_small, bed2_small, 
                   col = cols_vec_small)


# Add legend
library(plotrix)
color.legend(0.6, #left
             -1, #bottom
             0.9, #right
             -0.95, #top
             legend = c('50%', '100%'),
             align = 'rb',
             #rev(ddf$VAL),
             rect.col = (rbPal(20)),
             gradient="x",
             cex = 0.7,
             )
text(x = 0.75, y = -0.92,
     labels = 'similarity',
     cex = 0.7)





# Other computations ------------------------------------------------------
## Further operations not required for plotting, but that use the same inputs as for circos plot

## Get prop of each chr covered by a block ---------------------------------
#less "SaFo_self_masked_mindots30_topn2_blocks_interchr_blocks" | cut -f1,4,5 | sort -k1,1 -k2,2 -n > "SaFo_self_masked_mindots30_topn2_blocks_interchr_blocks.sorted.bed"
#bedtools merge -i "SaFo_self_masked_mindots30_topn2_blocks_interchr_blocks.sorted.bed" > "SaFo_self_masked_mindots30_topn2_blocks_interchr_blocks.sorted.merged.bed"
blocks.sorted.merged <- read.delim("~/SaFo_paper/SaFo_self_masked_mindots30_topn2_blocks_interchr_blocks.sorted.merged.bed", header=FALSE,
                                   col.names = c('CHROM', 'START', 'STOP'))
blocks.sorted.merged$BP <- blocks.sorted.merged$STOP - blocks.sorted.merged$START

blocks.sorted.merged <- merge(x = blocks.sorted.merged, y = rchr_lengths, 
                              by.x = 'CHROM', by.y = 'chrNum', sort = FALSE, all = TRUE)
blocks.sorted.merged_bp <- blocks.sorted.merged %>% group_by(CHROM, size) %>% summarise(BP = sum(BP))

blocks.sorted.merged_bp$prop <- blocks.sorted.merged_bp$BP/blocks.sorted.merged_bp$size


write.table(file = "~/SaFo_paper/SaFo_self_masked_mindots30_topn2_blocks_interchr_blocks_chr_props.txt",
            blocks.sorted.merged_bp,
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")



## High similarity regions (>90%) ---------------------------------------
high_sim <- subset(identity_win,IDY_W > 90) 


## Putatively collapsed regions -----------------------------------------
intersect_depth_syn <- read.delim("~/SaFo_paper/intersect_depth_synblocks_1Mb.txt", 
                                  header=FALSE, 
                                  col.names = c('CHROM', 'START', 'STOP', 'MEAN_DEPTH', 
                                                'BLOCK_CHROM', 'BLOCK_START', 'BLOCK_STOP', 'OVERLAP'))
intersect_depth_syn$high_depth <- ifelse(intersect_depth_syn$MEAN_DEPTH > mean(intersect_depth_syn$MEAN_DEPTH) + 2*(sd(intersect_depth_syn$MEAN_DEPTH)),
                         yes = 'yes', 
                         no = 'no')

intersect_depth_syn$homology <- ifelse(intersect_depth_syn$OVERLAP > 100,
                                                yes = 'yes',
                                                no = 'no')

intersect_depth_syn$collapsed <- ifelse(intersect_depth_syn$high_depth == 'yes' & intersect_depth_syn$homology == 'no',
                                        yes = 'yes',
                                        no = 'no')
collapsed <- subset(intersect_depth_syn, collapsed == 'yes')


write.table(file = "~/SaFo_paper/intersect_depth_synblocks_1Mb_putative_collapsed.txt",
            collapsed,
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


## Compute number of bp in putatively collapsed regions
sum(collapsed$STOP - collapsed$START) / TOT_CHR_SIZE





