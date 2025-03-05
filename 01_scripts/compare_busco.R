

safo <- read.delim("species_comparison/busco/BUSCO_busco/safo.fna/run_actinopterygii_odb10/full_table.tsv", header=FALSE, comment.char="#")[, c(1,2)]
safo <- safo[!duplicated(safo), ]
safo$sp <- 'safo'

sasa <- read.delim("species_comparison/busco/BUSCO_busco/sasaICSASG.fna/run_actinopterygii_odb10/full_table.tsv", header=FALSE, comment.char="#")[, c(1,2)]
sasa <- sasa[!duplicated(sasa), ]
sasa$sp <- 'sasa'

omyk <- read.delim("species_comparison/busco/BUSCO_busco/omyk.fna/run_actinopterygii_odb10/full_table.tsv", header=FALSE, comment.char="#")[, c(1,2)]
omyk <- omyk[!duplicated(omyk), ]
omyk$sp <- 'omyk'

sana <- read.delim("species_comparison/busco/BUSCO_busco/sana.fna/run_actinopterygii_odb10/full_table.tsv", header=FALSE, comment.char="#")[, c(1,2)]
sana <- sana[!duplicated(sana), ]
sana$sp <- 'sana'

cocl <- read.delim("species_comparison/busco/BUSCO_busco/cocl.fna/run_actinopterygii_odb10/full_table.tsv", header=FALSE, comment.char="#")[, c(1,2)]
cocl <- cocl[!duplicated(cocl), ]
cocl$sp <- 'cocl'

satr <- read.delim("species_comparison/busco/BUSCO_busco/satr.fna/run_actinopterygii_odb10/full_table.tsv", header=FALSE, comment.char="#")[, c(1,2)]
satr <- satr[!duplicated(satr), ]
satr$sp <- 'satr'

sama <- read.delim("species_comparison/busco/BUSCO_busco/sama.fna/run_actinopterygii_odb10/full_table.tsv", header=FALSE, comment.char="#")[, c(1,2)]
sama <- sama[!duplicated(sama), ]
sama$sp <- 'sama'



combined <- dplyr::bind_rows(list(safo, sana, sama, omyk, satr, cocl, sasa))

combined_wide <- tidyr::pivot_wider(combined, 
                                    names_from = sp,
                                    values_from = V2)

# Count number of Complete or Duplicated for each gene across species
combined_wide$Complete <- NA
combined_wide$Duplicated <- NA
combined_wide$Missing <- NA
for (i in 1:nrow(combined_wide)){
  combined_wide$Complete[i] <- length(grep(combined_wide[i, 2:8], pattern = 'Complete'))
  combined_wide$Duplicated[i] <- length(grep(combined_wide[i, 2:8], pattern = 'Duplicated'))
  combined_wide$Missing[i] <- length(grep(combined_wide[i, 2:8], pattern = 'Missing'))
  #combined_wide$Duplicated[i] <- length(grep(combined_wide[i, 3:9], pattern = 'Duplicated'))
  #combined_wide$Complete <- length(paste(combined_wide[i, 3:ncol(combined_wide)]) == 'Complete')
}


#combined_wide <- merge(combined_wide, links_to_ODB10[, 1:2], by = 'V1')

# Check complete across all sp
all_complete <- subset(combined_wide, Complete == 7)
all_duplicated <- subset(combined_wide, Duplicated == 7)
all_missing <- subset(combined_wide, Missing == 7)

# Check missing and duplicated in safo 
missing_safo <- subset(combined_wide, safo == 'Missing')
duplicated_safo <- subset(combined_wide, safo == 'Duplicated')


# Export 
write.table(combined_wide, file = 'species_comparison/busco/combined_busco.txt', sep = "\t", quote = FALSE, row.names = FALSE)






# Check duplicated across at least 3 sp
duplicated_over2 <- subset(combined_wide, Duplicated >= 3)



combined %>% group_by(V1) %>% 
  
  
  missing <- subset(full_table_busco_format, V2 == 'Missing')

merge(missing, links_to_ODB10, by = 'V1')