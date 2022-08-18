# Aim: To correlate the RSK4 isoform expression with overall survival of cancer patients

# Load packages 
library(survminer)
library(survival) # survival analysis 
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)

# Import data 
setwd("/Users/chensisi/Documents/RNAseq/")
load("merge_combine.RData")

# Univariate survival analysis by TPM

isoform_tpm <- merge_combine %>%
  select("ENST00000620340.4","ENST00000262752.4")


#### Survival analysis (in general) ####
# fit a univariate model for each isoform fraction, and find those most significantly associated to patient outcome
results.univariate <-array(NA, c(ncol(isoform_tpm),4)) # isoform_tpm / isoform_ratio 
colnames(results.univariate)<-c("HR","LCI","UCI","PVAL")
rownames(results.univariate) <-colnames(isoform_tpm) # isoform_tpm / isoform_ratio 
results.univariate <-as.data.frame(results.univariate)

for(i in 1:ncol(isoform_tpm)) # isoform_tpm / isoform_ratio
{
  coxphmodel <- coxph(merge2.os ~ as.numeric(isoform_tpm[,i])) # isoform_tpm / isoform_ratio
  results.univariate$HR[i] <- summary(coxphmodel)$coef[1,2]
  results.univariate$LCI[i] <- summary(coxphmodel)$conf.int[1,3]
  results.univariate$UCI[i] <- summary(coxphmodel)$conf.int[1,4]
  results.univariate$PVAL[i] <- summary(coxphmodel)$coef[1,5]
}

# Calculate P-adjusted value (FDR) and re-order the data.frame by FDR value 
results.univariate$FDR<- p.adjust(results.univariate$PVAL,method="fdr")
results.univariate<-results.univariate[order(results.univariate$FDR, decreasing=F),]
results.univariate



#### Univariate Survival analysis (by cancer type) ####
cancer_types <- unique(merge_combine$Sample_type) 
my_list <- list() # create a list to store results of survival analysis 

univariate.results <- data.frame()

for (x in cancer_types){
  cancer <- merge_combine %>% 
    filter(Sample_type == x)

  # Survival analysis 
  cancer.os <- Surv(cancer$OS.time, cancer$OS)
  cancer.IF <- cancer %>%
    select("ENST00000620340.4","ENST00000262752.4") # "ratio"
  
  results <- array(NA, c(ncol(cancer.IF),4))
  colnames(results)<-c("HR","LCI","UCI","PVAL")
  rownames(results)<-colnames(cancer.IF)
  results <- as.data.frame(results)
  
  for(i in 1:ncol(cancer.IF)) 
  {
    coxphmodel <- coxph(cancer.os ~ as.numeric(cancer.IF[,i]))
    
    results$HR[i] <- summary(coxphmodel)$coef[1,2]
    
    results$LCI[i] <- summary(coxphmodel)$conf.int[1,3]
    
    results$UCI[i] <- summary(coxphmodel)$conf.int[1,4]
    
    results$PVAL[i] <- summary(coxphmodel)$coef[1,5]
  }
  
  results$FDR <- p.adjust(results$PVAL, method ="fdr")
  
  results$cancer_type <- x

  my_list[[x]] <- results
  
  print(paste0("Survival analysis for ",x, " has finished!"))

}

univariate.results_tpm <- data.frame()
for (i in 1:length(my_list)) {
  survival <- as.data.frame(my_list[[i]])
  survival <- cbind(rownames(survival), survival)
  rownames(survival) <- NULL
  univariate.results <- rbind(univariate.results, survival)
}

univariate.results <- univariate.results[order(univariate.results$FDR, decreasing=F),]
univariate.results <- univariate.results %>% 
  select("rownames(survival)","HR","LCI","UCI","PVAL","FDR","cancer_type" )

univariate.results[1:10,]

# Export results 
write.table(univariate.results, "univariate.results.csv", sep =",", col.names = TRUE, row.names = FALSE)



