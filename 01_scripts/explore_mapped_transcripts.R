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


# Plot distrib
ggplot(data = mapped_uniques_80) + 
  facet_wrap(.~species) + 
  geom_histogram(aes(x = identity), binwidth = 1) + 
  lims(x = c(60, 102))

# Plot copy number within asms
mapped_barplot <- 
ggplot(data = mapped_uniques_80, aes(x = factor(species, level = c('safo', 'sana', 'sama', 'omyk', 'satr', 'sasaICSASG', 'cocl')))) + 
  geom_bar(aes(fill = factor(cp_number))) +
  scale_fill_viridis_d(option = 'H') + 
  labs(fill='Transcript copy number', x = 'Species', y = 'Transcript count') +
  scale_x_discrete(labels = c(
    "safo" = "Salvelinus fontinalis", 
    "sana" = "Salvelinus namaycush", 
    "sama" = "Salvelinus sp.",
    "omyk" = "Oncorhynchus mykiss",
    "satr" = "Salmo trutta",
    "sasaICSASG" = "Salmo salar",
    "cocl" = "Coregonus clupeaformis")) +
  scale_y_continuous(labels = scales::number_format(big.mark = ",")) +
  theme(
    axis.text.x = element_text(face = 'italic', angle = 45, size = 9, hjust = 1),
    axis.text.y = element_text(size = 9, hjust = 1),
    legend.title = element_text(size = 9),
    panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.2),
    panel.background = element_blank()
  )

ggsave(filename = 'species_comparison/transcriptome/mapped_transcripts_over80idy_barplot.pdf',
       mapped_barplot, width = 10, height = 8)

#missing_safo_transcripts <- subset(mapped_uniques_wide_80, safo == 0)
#duplicated_safo <- subset(mapped_uniques_wide, safo > 1)



# Export
write.table(mapped_uniques_wide_80, file = 'species_comparison/transcriptome/mapped_transcripts_over80idy.txt', sep = "\t", quote = FALSE, row.names = FALSE)



