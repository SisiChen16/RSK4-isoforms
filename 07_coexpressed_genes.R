library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(EnvStats)

##### Import Data ####
setwd("/Users/chensisi/Desktop/MRes_Project1_RSK4/Codes/Data/")

# Transcript data import 
transcript <- read.table("TCGA_GTEx_RSK4_tpm.txt", header = TRUE,  row.names=1)
transcript1 <- as.data.frame(t(transcript)) # Transpose the data frame 
transcript1 <- transcript1[,c("ENST00000620340.4","ENST00000262752.4")]
colnames(transcript1) <- c("Iso1","Iso2")

LGG <- read.table("LGGexpression.tpm.txt", header = TRUE, row.names = 1, sep = "\t")
STAD <- read.table("STADexpression.tpm.txt", header = TRUE, row.names = 1, sep = "\t")
KIRC <- read.table("KIRCexpression.tpm.txt", header = TRUE, row.names = 1, sep = "\t")
CESC <- read.table("CESCexpression.tpm.txt", header = TRUE, row.names = 1, sep = "\t")
# Convert the dataframe into matrix 
LGG <- as.data.frame(LGG)
STAD <- as.data.frame(STAD)
KIRC <- as.data.frame(KIRC)
CESC <- as.data.frame(CESC)

# Metadata import 
metadata <- read.table("TCGA_GTEX_category.txt",  header = TRUE, row.names = 1, fill = TRUE,sep = "\t")
metadata <- metadata %>% separate(TCGA_GTEX_main_category, c("Project", "Sample_type"), sep = "\\s", extra = "merge")
rownames(metadata)<-gsub(rownames(metadata), pattern="-", replace=".")
metadata[1:3,]

## Merge multiple data frames by matching sample IDs 
merge1 <- merge(transcript1, metadata,  by = 0)
rownames(merge1) <- merge1[,1]
merge1 <- merge1[,-1]

brain <- merge1 %>% filter(Project == "TCGA", Sample_type == "Brain Lower Grade Glioma")
stomach <- merge1 %>% filter(Project == "TCGA", Sample_type == "Stomach Adenocarcinoma")
kidney <- merge1 %>% filter(Project == "TCGA", Sample_type == "Kidney Clear Cell Carcinoma")
cervical <- merge1 %>% filter(Project == "TCGA", Sample_type == "Cervical & Endocervical Cancer")


## Prepare the dataframes
  # Brain LGG 
brain2 <- as.data.frame(t(brain))
brain2 <- brain2[,which(is.element(colnames(brain2),colnames(LGG)))]
brain2 <- brain2[,order(colnames(brain2))]
LGG2 <- LGG[,order(colnames(LGG))]
head(colnames(brain2))
head(colnames(LGG2))

  # Stomach STAD 
stomach2 <- as.data.frame(t(stomach))
stomach2 <- stomach2[,which(is.element(colnames(stomach2),colnames(STAD)))]
stomach2 <- stomach2[,order(colnames(stomach2))]
STAD2 <- STAD[,order(colnames(STAD))]
head(colnames(stomach2))
head(colnames(STAD2))

  # Kidney KIRC 
kidney2 <- as.data.frame(t(kidney))
kidney2 <- kidney2[,which(is.element(colnames(kidney2),colnames(KIRC)))]
KIRC2 <- KIRC[,which(is.element(colnames(KIRC),colnames(kidney2)))]
kidney2 <- kidney2[,order(colnames(kidney2))]
KIRC2 <- KIRC[,order(colnames(KIRC2))]
head(colnames(kidney2))
head(colnames(KIRC2))

  # Cervical CESC 
cervical2 <- as.data.frame(t(cervical))
cervical2 <- cervical2[,which(is.element(colnames(cervical2),colnames(CESC)))]
CESC2 <- CESC[,which(is.element(colnames(CESC),colnames(cervical2)))]
cervical2 <- cervical2[,order(colnames(cervical2))]
CESC2 <- CESC[,order(colnames(CESC2))]
head(colnames(cervical2))
head(colnames(CESC2))

# Correlation tests 

  # Brain - LGG 
LGG.results <-array(NA,c(nrow(LGG2),2))
rownames(LGG.results)<-rownames(LGG2)
colnames(LGG.results)<- c("rho","Pval")

for (i in 1:nrow(LGG2)){
  LGG.results[i,1] <- cor.test(as.numeric(brain2["Iso1",]), as.numeric(LGG2[i,]), 
                           method  = "spearman")$est 
  
  LGG.results[i,2] <- cor.test(as.numeric(brain2["Iso1",]), as.numeric(LGG2[i,]), 
                           method  = "spearman")$p.value}
LGG.results <- as.data.frame(LGG.results)
LGG.results$p.adj <- p.adjust(LGG.results$Pval, method = "fdr")

  # Stomach - STAD 
STAD.results<-array(NA,c(nrow(STAD2),2))
rownames(STAD.results)<-rownames(STAD2)
colnames(STAD.results)<- c("rho","Pval")

for (i in 1:nrow(STAD2)){
  STAD.results[i,1] <- cor.test(as.numeric(stomach2["Iso1",]), as.numeric(STAD2[i,]), 
                               method  = "spearman")$est 
  
  STAD.results[i,2] <- cor.test(as.numeric(stomach2["Iso1",]), as.numeric(STAD2[i,]), 
                               method  = "spearman")$p.value}
STAD.results <- as.data.frame(STAD.results)
STAD.results$p.adj <- p.adjust(STAD.results$Pval, method = "fdr")

  # Kidney - KIRC 
KIRC.results<-array(NA,c(nrow(KIRC2),2))
rownames(KIRC.results)<-rownames(KIRC2)
colnames(KIRC.results)<- c("rho","Pval")

for (i in 1:nrow(KIRC2)){
  KIRC.results[i,1] <- cor.test(as.numeric(kidney2["Iso1",]), as.numeric(KIRC2[i,]), 
                                method  = "spearman")$est 
  
  KIRC.results[i,2] <- cor.test(as.numeric(kidney2["Iso1",]), as.numeric(KIRC2[i,]), 
                                method  = "spearman")$p.value}
KIRC.results <- as.data.frame(KIRC.results)
KIRC.results$p.adj <- p.adjust(KIRC.results$Pval, method = "fdr")

  # Cervical - CESC
CESC.results1<-array(NA,c(nrow(CESC2),2))
rownames(CESC.results1)<-rownames(CESC2)
colnames(CESC.results1)<- c("rho","Pval")

for (i in 1:nrow(CESC2)){
  CESC.results1[i,1] <- cor.test(as.numeric(cervical2["Iso1",]), as.numeric(CESC2[i,]), 
                                method  = "spearman")$est 
  
  CESC.results1[i,2] <- cor.test(as.numeric(cervical2["Iso1",]), as.numeric(CESC2[i,]), 
                                method  = "spearman")$p.value}
CESC.results1 <- as.data.frame(CESC.results1)
CESC.results1$p.adj <- p.adjust(CESC.results1$Pval, method = "fdr")

# Data processing 
LGG.results <- LGG.results %>% arrange(desc(rho)) %>% filter(!is.na(rho))
STAD.results <- STAD.results %>% arrange(desc(rho)) %>% filter(!is.na(rho))
KIRC.results <- KIRC.results %>% arrange(desc(rho)) %>% filter(!is.na(rho))
CESC.results <- CESC.results %>% arrange(desc(rho)) %>% filter(!is.na(rho))

#write.table(LGG.results, "LGG.results.full.csv", sep = ",")
#write.table(STAD.results,"STAD.results.full.csv",sep = ",")
#write.table(KIRC.results, "KIRC.results.full.csv",sep = ",")
#write.table(CESC.results, "CESC.results.full.csv",sep = ",")


KIRC.results1 <- KIRC.results1 %>% arrange(desc(rho)) %>% filter(!is.na(rho))
CESC.results1 <- CESC.results1 %>% arrange(desc(rho)) %>% filter(!is.na(rho))
#write.table(KIRC.results1, "KIRC.results1.full.csv",sep = ",")
#write.table(CESC.results1, "CESC.results1.full.csv",sep = ",")

# Filter out significant results 
LGG.sig <- LGG.results %>% filter(rho > 0.55, p.adj < 0.05) %>% arrange(desc(rho))
LGG.coexpress <- rownames(LGG.sig)[1:101]
LGG.sig <- LGG.results %>% filter(rho < -0.4, p.adj < 0.05) %>% arrange(rho)
LGG.repress <- rownames(LGG.sig)[1:100]
#write.table(LGG.coexpress, "LGG.coexpress.100.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\n")

STAD.sig <- STAD.results %>% filter(rho > 0.55, p.adj < 0.05) %>% arrange(desc(rho))
STAD.coexpress <- rownames(STAD.sig)[1:101]
STAD.sig <- STAD.results %>% filter(rho < -0.35, p.adj < 0.05) %>%  arrange(rho) # -0.38 
STAD.repress <- rownames(STAD.sig)[1:100]
#write.table(STAD.repress, "STAD.repress.100.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\n")

KIRC.results1 %>% filter(abs(rho) > 0.3 , p.adj < 0.05)

CESC.sig <- CESC.results1 %>% filter(abs(rho) > 0.3 , p.adj < 0.05 )
CESC.coexpress <- rownames(CESC.sig)
#write.table(CESC.coexpress, "CESC.coexpress.66.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\n")
