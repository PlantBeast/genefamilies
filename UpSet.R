library("UpSetR")

#Set path
setwd("/Users/sstrickler/Desktop/Basil/Orthogroups/")

orthogroups_df<- read.table("Orthogroups.GeneCount.tsv",  header=T, stringsAsFactors = F)

#All species
selected_species <- colnames(orthogroups_df)[2:(ncol(orthogroups_df) -1)] 
selected_species
ncol(orthogroups_df)
orthogroups_df[orthogroups_df > 0] <- 1
upset(orthogroups_df, nsets = ncol(orthogroups_df), sets = rev(selected_species), 
      keep.order = T, order.by = "freq", number.angles = 30, empty.intersections = "on")
