library(ggplot2)
library(dplyr)

# Import
## Ssal_v3.1 transcripts
# Transcripts
Ssal_v3.1_rna <- read.delim("species_comparison/transcriptome/GCF_905237065.1_Ssal_v3.1_rna.txt", header=FALSE,
                            col.names = c('name', 'full_name'))
# Mapped transcripts for each asm
for (i in c('safo', 'sana', 'sama', 'omyk', 'satr', 'sasaICSASG', 'cocl')){
  assign(x = eval(i), 
         value = read.delim(
           paste0("species_comparison/transcriptome/", i, '/Ssal_v3.1_on_', i,'_gmap_mRNA.txt'), header=FALSE, 
           col.names = c('chr', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'ID', 'name', 'parent', 'coverage', 'identity', 'matches', 'mismatches', 'indels', 'unknowns')
           )
  )
}

# Filter mapped transcripts on % identity and get copy number 
dfs <- list(safo, sana, sama, omyk, satr, sasaICSASG, cocl)
df_names <- c('safo', 'sana', 'sama', 'omyk', 'satr', 'sasaICSASG', 'cocl')

for (i in 1:length(dfs)){
  
  # Filter on identity or other criteria if required
  # dfs[[i]] <- subset(dfs[[i]], identity >= 80)
  ### transcripts with identity <80 will be treated as unmapped/missing
  
  # Add copy number 
  counts_df <- dfs[[i]] %>% count(name)
  df_merged <- merge(x = dfs[[i]], y = counts_df, by = 'name')
  
  # Add present-absent info  
  df_merged <- merge(x = df_merged, y = Ssal_v3.1_rna, by = 'name', all.y = TRUE) # faire avant ajout variable species
  df_merged$status <- ifelse(test = is.na(df_merged$identity), yes = 'absent', no = 'present')
  
  # Add copy number = 0 if absent
  df_merged$cp_number <- ifelse(test = df_merged$status == 'absent',
                                  yes = 0, no = df_merged$n)
  # Add variable for species
  df_merged$species <- df_names[i]
  
  
  dfs[[i]] <- df_merged
} 

# Combine all species into a single df
all_mapped <- do.call(rbind, dfs)

# Keep only one occurence of each transcript
mapped_uniques <- all_mapped[!duplicated(all_mapped[, c('name', 'species')]), ]

# Plot --------------------------------------------------------------------
ggplot(data = mapped_uniques) + 
  facet_wrap(.~species) + 
  geom_histogram(aes(x = identity), binwidth = 1) + 
  lims(x = c(60, 102))

ggplot(data = mapped_uniques) + 
  geom_bar(aes(x = species, fill = factor(cp_number)))



# Summarize
summary_mRNA <- function(sp, x) {
  unmapped <- length(unique(x[x$status == 'absent', 'name']))
  mapped <- length(unique(x[x$status == 'present', 'name']))
  duplicated <- length(unique(x[x$cp_number > 1, 'name']))
  
  print(paste(sp, ': unmapped =', unmapped, '; mapped =', mapped, '; duplicated =', duplicated))
} 

summary_idy <- function(sp, x) {
  print(paste(sp, ': min =', min(x$identity, na.rm = TRUE), '; max =', max(x$identity, na.rm = TRUE), 
              '; mean =', mean(x$identity, na.rm = TRUE), '; mediam = ', median(x$identity, na.rm = TRUE) ))
} 


for (sp in c('safo', 'sana', 'sama', 'omyk', 'satr', 'sasaICSASG', 'cocl') ){
  summary_mRNA(sp, subset(mapped_uniques, species == sp))
  summary_idy(sp, subset(mapped_uniques, species == sp))
}


# Get copy number per transcript per species
mapped_uniques_wide <- tidyr::pivot_wider(mapped_uniques[!duplicated(mapped_uniques[, c('name', 'species')]), c('name', 'species',  'cp_number')], 
                                      names_from = species,
                                      values_from = cp_number)



# With 80% cutoff ---------------------------------------------------------

# Add infos
dfs_80 <- list(safo, sana, sama, omyk, satr, sasaICSASG, cocl)
#df_names <- c('safo', 'sana', 'sama', 'omyk', 'satr', 'sasaICSASG', 'cocl')

for (i in 1:length(dfs_80)){
  
  # Filter on identity or other criteria if required
  dfs_80[[i]] <- subset(dfs_80[[i]], identity >= 80)
  ### transcripts with identity <80 will be treated as unmapped/missing
  
  # Add copy number 
  counts_df <- dfs_80[[i]] %>% count(name)
  df_merged <- merge(x = dfs_80[[i]], y = counts_df, by = 'name')
  
  # Add present-absent info  
  df_merged <- merge(x = df_merged, y = Ssal_v3.1_rna, by = 'name', all.y = TRUE) # faire avant ajout variable species
  df_merged$status <- ifelse(test = is.na(df_merged$identity), yes = 'absent', no = 'present')
  
  # Add cp number = 0 if absent
  df_merged$cp_number <- ifelse(test = df_merged$status == 'absent',
                                yes = 0, no = df_merged$n)
  # Add variable for species
  df_merged$species <- df_names[i]
  
  
  dfs_80[[i]] <- df_merged
} 

# Combine all species into a single df
all_mapped_80 <- do.call(rbind, dfs_80)

# Keep only one occurence of each transcript
mapped_uniques_80 <- all_mapped_80[!duplicated(all_mapped_80[, c('name', 'species')]), ]

# Get copy number per transcript per species
mapped_uniques_wide_80 <- tidyr::pivot_wider(mapped_uniques_80[!duplicated(mapped_uniques_80[, c('name', 'species')]), c('name', 'species',  'cp_number')], 
                                          names_from = species,
                                          values_from = cp_number)

# Get distance between species
out_dist <- expand.grid('species1' = df_names, 'species2' = df_names)
for (i in 1:nrow(out_dist)){
  out_dist$dist[i] <- sum(abs(mapped_uniques_wide_80[[out_dist$species1[i]]] - mapped_uniques_wide_80[[out_dist$species2[i]]]))/nrow(mapped_uniques_wide_80)
}

tidyr::pivot_wider(out_dist, 
                   names_from = c('species1'),
                   values_from = 'dist')



# Plot distrib
ggplot(data = mapped_uniques_80) + 
  facet_wrap(.~species) + 
  geom_histogram(aes(x = identity), binwidth = 1) + 
  lims(x = c(60, 102))

# Plot copy number within asms
ggplot(data = mapped_uniques_80) + 
  geom_bar(aes(x = species, fill = factor(cp_number))) +
  scale_fill_viridis_d(option = 'D')



mapped_dist <- dist(mapped_uniques_wide_80[, 2:8])


missing_transcripts_safo <- subset(mapped_uniques_wide_80, safo == 0)
duplicated_transcripts_safo <- subset(mapped_uniques_wide_80, safo > 1)



for (i in 1:nrow(mapped_uniques_wide_80)){
  mapped_uniques_wide_80$single_copy[i] <- sum(mapped_uniques_wide_80[i, 2:8] == 1)
  mapped_uniques_wide_80$duplicated[i] <- sum(mapped_uniques_wide_80[i, 2:8] > 1)
  mapped_uniques_wide_80$missing[i] <- sum(mapped_uniques_wide_80[i, 2:8] == 0)
}

missing_safo_transcripts <- subset(mapped_uniques_wide_80, safo == 0)
duplicated_safo <- subset(mapped_uniques_wide, safo > 1)



# Export
write.table(mapped_uniques_wide_80, file = 'species_comparison/transcriptome/mapped_transcripts_over80idy.txt', sep = "\t", quote = FALSE, row.names = FALSE)
















all_mapped_wide <- tidyr::pivot_wider(all_mapped[!duplicated(all_mapped[, c('name', 'species')]), c('name', 'species',  'cp_number')], 
                                    names_from = species,
                                    values_from = cp_number)



out_df <- data_frame(name = unique(Ssal_v3.1_rna$name))

for (sp in c('safo', 'sana', 'sama', 'omyk', 'satr', 'sasaICSASG', 'cocl')){
  for (i in unique(Ssal_v3.1_rna$name)){
    df <- subset(all_mapped, name == i & species == sp)
    if (all(df$cp_number == 1)) {
      print(paste0('single copy (', df$identity, ')'))
    } else if (nrow(df) > 1) {
      print(paste0('duplicated (', cat(df$identity, sep = ','), ')'))
    } else if (df$cp_number == 0){
      print('missing')
    }
    }
}

out_df <- data_frame(name = unique(Ssal_v3.1_rna$name))












for (sp in c('safo', 'sana', 'sama', 'omyk', 'satr', 'sasaICSASG', 'cocl')){
  
  #assign(x = paste0(sp, '_status'), value = vector(mode = 'character', length = nrow(Ssal_v3.1_rna)))
  status_vec <- vector(mode = 'character', length = nrow(Ssal_v3.1_rna))
  idy_vec <- vector(mode = 'character', length = nrow(Ssal_v3.1_rna))
  sp_df <- subset(all_mapped, species == sp)
  
  for (i in unique(Ssal_v3.1_rna$name)){
    
    df <- subset(sp_df, name == i)
    #print(cat(df$identity, sep = ','))
    status_vec[i] <- 
      (ifelse(test = all(df$cp_number == 1), yes = 'single copy', 
                            no = ifelse(test = all(df$cp_number > 1), yes = 'duplicated', no = 'missing')))
    idy_vec <- 
      (paste(df$identity))
  }
  assign(x = paste0('status_', sp), value = status_vec)
  assign(x = paste0('idy_', sp), value = idy_vec)
  }
















all_mapped_wider <- tidyr::pivot_wider(all_mapped[, c('name', 'species',  'cp_number', 'identity')], 
                                      names_from = species,
                                      values_from = cp_number)


tidyr::pivot_wider(all_mapped[, c('name', 'species',  'cp_number', 'identity')], 
                   names_from = species,
                   values_from = c(cp_number, identity)
                   )


mapped_uniques_wide$single_copy <- NA


safo_s <- subset(all_mapped, species == 'safo')

for (i in unique(safo_s$name)){
  if (safo_s$cp_number[safo_s$name == i] == 1){
    
    print(paste0('Single copy (', safo_s$identity[safo_s$name == i], ')'))
  } 
}

for (i in 1:nrow(all_mapped)){
  if (all_mapped$cp_number[i] == 1){
    print(paste0('single copy (', all_mapped$identity, ')'))
  }
  
}



# Poubelle ----------------------------------------------------------------


# Transcripts
Ssal_v3.1_rna <- read.delim("species_comparison/transcriptome/GCF_905237065.1_Ssal_v3.1_rna.txt", header=FALSE,
                                            col.names = c('name', 'full_name'))

length(setdiff(Ssal_v3.1_rna$name, sasaICSASG$name))
# 185 unmapped transcripts 
# 112712 mapped transcript names

length(intersect(Ssal_v3.1_rna$name, dfs[[5]]$name))
length(setdiff(Ssal_v3.1_rna$name, dfs[[5]]$name))

# add copy number info
satr_counts <- dfs[[5]] %>% count(name)
satr_merged <- merge(x = dfs[[5]], y = satr_counts, by = 'name')

# add present-absent info  
satr_merged <- merge(x = satr_merged, y = Ssal_v3.1_rna, by = 'name', all.y = TRUE) # faire avant ajout variable species
satr_merged$status <- ifelse(test = is.na(satr_merged$identity), yes = 'absent', no = 'present')
## add cp number = 0 if absent
satr_merged$cp_number <- ifelse(test = satr_merged$status == 'absent',
                                yes = 0, no = satr_merged$n)




dfs[[5]][ dfs[[5]]$name %in% dfs[[5]]$name[duplicated(dfs[[5]]$name)], ] # duplicated
dfs[[5]][ ! dfs[[5]]$name %in% dfs[[5]]$name[duplicated(dfs[[5]]$name)], ] # single copy

dfs[[5]][ dfs[[5]]$name %in% setdiff(Ssal_v3.1_rna$ID, dfs[[5]]$name), ]




library(ggplot2)


for (i in dfs){
  print(median(dfs[['identity']]))
}
