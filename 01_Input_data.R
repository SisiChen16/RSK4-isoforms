##　Project Aim: Differential expression analysis of RSK4 isoforms in tumour and normal tissue samples using TCGA and GTEx data


## Load Packages 
library(dplyr)
library(tidyverse)


##　Import data 
setwd("/Users/chensisi/Documents/RNAseq/")

# Transcript data import 
transcript <- read.table("TCGA_GTEx_RSK4_tpm.txt", header = TRUE,  row.names=1)
transcript1 <- as.data.frame(t(transcript)) # Transpose the data frame 
transcript1[1:3,]

# Metadata import 
metadata <- read.table("TCGA_GTEX_category.txt",  header = TRUE, row.names = 1, fill = TRUE,sep = "\t")

metadata <- metadata %>% 
  separate(TCGA_GTEX_main_category, c("Project", "Sample_type"), sep = "\\s", extra = "merge")

rownames(metadata)<-gsub(rownames(metadata), pattern="-", replace=".")
metadata[1:3,]

# Phenotype data import 
phenotype <- read.table("TcgaTargetGTEX_phenotype.txt", sep="\t", header=T, row.names=1) 
rownames(phenotype)<-gsub(rownames(phenotype), pattern="-", replace=".")
phenotype[1:3,]


# Survival data import 
survival <- read.table("TCGA_survival_data", sep="\t", header=T, row.names=1) 
rownames(survival)<-gsub(rownames(survival), pattern="-", replace=".")
survival[1:3,]

# PANCAN phenotype data 
pancan <- read.table("PANCAN.phenotype.txt", sep="\t", header=T, row.names=1)
rownames(pancan)<-gsub(rownames(pancan), pattern="-", replace=".")
colnames(pancan)
pancan <- pancan %>%  # remove the OS and OS.time column which duplicate with survival data
  select(-c(OS, OS.time))
pancan[1:3,]


## Merge multiple data frames by matching samples IDs 
merge1 <- merge(transcript1, metadata,  by = 0) # transcript data merge with metadata 
rownames(merge1) <- merge1[,1]
merge1 <- merge1[,-1]

merge1 <- merge(merge1, phenotype, by = 0) # merge with phenotype data
rownames(merge1) <- merge1[,1]
merge1 <- merge1[,-1]

merge2 <- merge(merge1, survival, by = 0) # merge with survival data
rownames(merge2) <- merge2[,1]
merge2 <- merge2[,-1]

merge3 <- merge(merge2, pancan, by = 0) # merge with pancancer phenotype/clinical data 
rownames(merge3) <- merge3[,1]
merge3 <- merge3[,-1]


##　Calculate the ratio of isoform2/isoform1 for later comparison 

merge_combine <- merge3 %>%
	mutate(iso1_denorm = 2^ENST00000620340.4,
		   iso2_denorm = 2^ENST00000262752.4) %>%
	mutate(ratio = iso2_denorm/iso1_denorm)

##　Split the data frame by project (TCGA or GTEX) & group them by sample types for later comparison

TCGA_by_sample <- merge_combine %>% 
  filter(Project == "TCGA") %>% 
  group_by(Sample_type) 

GTEX_by_sample <- merge_combine %>% 
  filter(Project == "GTEX") %>% 
  group_by(Sample_type)


## Save the data 
save.image(file = "RSK4_data.RData")
save(merge_combine, file = "merge_combine.RData") # save the data frame and reload it later


