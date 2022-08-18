# Aim: To test whether RSK4 isoform expression acts as an independnet prognostic predictor of patients' survival (by taking into consideration of other potential covariates)


# Load packages 
library(survminer)
library(survival) # survival analysis 
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)

#### Multivariate Survival analysis (by cancer type) with age + gender as covariates ####
cancer_types <- unique(merge_combine$Sample_type) 
my_list <- list() # create a list to store results of survival analysis 

coxph_list <- list()
for (x in cancer_types){
  cancer <- merge_combine %>% 
    filter(Sample_type == x)
  
  cancer.os <- Surv(cancer$OS.time, cancer$OS)
  
  # Potential factors
  cancer.pheno <- merge3 %>% 
    filter(Sample_type == x)

    # Pathological/clinical stage 
  cancer.pheno[cancer.pheno == ""] <- NA  # fill empty cell with NA
  cancer.pheno <- cancer.pheno %>%
    mutate(Merged_pathologic_clinical_stage = coalesce(ajcc_pathologic_tumor_stage, clinical_stage)) 

  pathologic_stage <- as.factor(cancer.pheno$Merged_pathologic_clinical_stage)
  x3<- grep("III",cancer.pheno$Merged_pathologic_clinical_stage)
  x4<- grep("IV",cancer.pheno$Merged_pathologic_clinical_stage)
  stage.high<-rep(0,nrow(cancer.pheno))
  stage.high[c(x3,x4)]<-1

  summary(coxph(cancer.os ~ stage.high))$coef

    # Gender 
  female <- grep("FEMALE", cancer.pheno$gender)
  gender <- rep(0,nrow(cancer.pheno))
  gender[female]<-1

  summary(coxph(cancer.os ~ gender)) 
  
    # Age 
  age <- cancer.pheno$age_at_initial_pathologic_diagnosis
  summary(coxph(cancer.os ~ age))


  # Cancer recurrence & distant metastasis 
  table(cancer.pheno$new_tumor_event_type)
  re <- grep("Recurrence",cancer.pheno$new_tumor_event_type)
  me <- grep("Metastasis|Metastatic",cancer.pheno$new_tumor_event_type)
  recur_metastasis <- rep(0,nrow(cancer.pheno))
  recur_metastasis[c(re,me)]<- 1

  summary(coxph(cancer.os ~ recur_metastasis))
  
  # Survival analysis 
  cancer.IF <- cancer %>%
    select("ENST00000620340.4","ENST00000262752.4") #"ENST00000620340.4","ENST00000262752.4","IF1","IF2"  "ratio"
  
  results <- array(NA, c(ncol(cancer.IF),4))
  colnames(results)<-c("HR","LCI","UCI","PVAL")
  rownames(results)<-colnames(cancer.IF)
  results <- as.data.frame(results)
  
  for(i in 1:ncol(cancer.IF)) 
  {
    coxphmodel <- coxph(cancer.os ~ cancer.IF[,i] + age + gender )
    
    results$HR[i] <- summary(coxphmodel)$coef[1,2]
    
    results$LCI[i] <- summary(coxphmodel)$conf.int[1,3]
    
    results$UCI[i] <- summary(coxphmodel)$conf.int[1,4]
    
    results$PVAL[i] <- summary(coxphmodel)$coef[1,5]
  }
  
  coxph_list[[x]] <- summary(coxphmodel)
  
  results$cancer_type <- x

  my_list[[x]] <- results
  
  print(paste0("Survival analysis for ",x, " has finished!"))

}


multivariate.results_tpm <- data.frame()
for (i in 1:length(my_list)) {
  survival <- as.data.frame(my_list[[i]])
  survival <- cbind(rownames(survival), survival)
  rownames(survival) <- NULL
  multivariate.results_tpm <- rbind(multivariate.results_tpm, survival)
}

multivariate.results_tpm$FDR<- p.adjust(multivariate.results_tpm$PVAL,method="fdr") 

multivariate.results_tpm <- multivariate.results_tpm[order(multivariate.results_tpm$PVAL, decreasing=F),]

multivariate.results_tpm[1:10,] 

multivariate.results_tpm <- multivariate.results_tpm %>% 
  select("rownames(survival)","HR","LCI","UCI","PVAL","FDR","cancer_type" )

# Export the multivariate survival results 
write.table(multivariate.results_tpm, "multivariate.results_by.tpm.csv", sep =",", col.names = TRUE, row.names = FALSE)


#### Take those cancer types which showed significant association between isoform expression and patients' survival
#### Re-do the multivariate survival analysis case by case by considering other potential variables (e.g. family history, headache, alcohol intake, etc)

### Multivariate survival analysis for cancer types of interest

# Import phenotype data for each cancer type
setwd("/Users/chensisi/Documents/RNAseq/3_Correlation_analysis/")

filenames <- c("TCGA.LGG.sampleMap_LGG_clinicalMatrix","TCGA.ACC.sampleMap-ACC_clinicalMatrix",
               "TCGA.STAD.sampleMap-STAD_clinicalMatrix","TCGA.KIRC.sampleMap_KIRC_clinicalMatrix",
               "TCGA.CESC.sampleMap_CESC_clinicalMatrix", "TCGA.READ.sampleMap_READ_clinicalMatrix") 

cancer.list <- c("LGG","ACC","STAD","KIRC","CESC","READ")

for (i in 1:length(filenames)){
  
  # Input cancer phenotype data for each cancer type 
  phenotype <- read.table(filenames[i], sep = '\t', header = TRUE, row.names = 1)
  rownames(phenotype)<-gsub(rownames(phenotype), pattern="-", replace=".")
  
  # Filter out the expression data for specific cancer type 
  expression <- merge_combine[which(is.element(rownames(merge_combine), rownames(phenotype))),]
  phenotype2 <- phenotype[which(is.element(rownames(phenotype),rownames(expression))),]
  
  # Re-order the dataframe
  expression <- as.data.frame(expression[order(rownames(expression)),])
  phenotype2 <-as.data.frame(phenotype2[order(rownames(phenotype2)),])
  
  # Assign the dataframe use different names
  assign(paste0(cancer.list[i],".exp"), data.frame(expression))
  assign(paste0(cancer.list[i],".pheno"), data.frame(phenotype2))
}


############################ Brain Lower Grade Glioma (LGG) ############################# 

load("merge_combine.RData")
LGG.rsk <- merge_combine %>% 
  filter(Sample_type == "Brain Lower Grade Glioma") %>% 
  select(ENST00000620340.4,ENST00000262752.4, ratio, OS, OS.time)

LGG.merge <- merge(LGG.rsk,LGG.pheno,by = 0)

# Create survival object 
LGG.os <- Surv(LGG.merge$OS.time, LGG.merge$OS)
# Covariate 1: family history of cancer 
family.history <- as.factor(LGG.merge$family_history_of_cancer)
table(family.history)
yes <- grep("YES",family.history)
family.history <- rep(0,nrow(LGG.merge))
family.history[yes] <- 1
summary(coxph(LGG.os ~ family.history)) #not significant (p = 0.283)

# Covariate 2: Gender 
gender <- as.factor(LGG.merge$gender)
table(gender)
summary(coxph(LGG.os ~ gender)) #Not significant (p = 0.538)

# Covariate 3: Age 
table(LGG.merge$age_at_initial_pathologic_diagnosis)
over60 <- which(LGG.merge$age_at_initial_pathologic_diagnosis > 60)
LGG.merge$age_over_60 <- 0 
LGG.merge$age_over_60[over60] <- 1
age <- LGG.merge$age_over_60
summary(coxph(LGG.os ~ age)) # Significant (p = 8.3e-14 ***)

# Covariate 4: Headache history
headache <- LGG.merge$headache_history
table(headache)
yes <- grep("YES",headache)
headache <- rep(0,nrow(LGG.merge))
headache[yes] <- 1
summary(coxph(LGG.os ~ headache)) # Not significant (p = 0.427)

# Covariate 5: Cancer recurrence 
recur <- LGG.merge$new_tumor_event_after_initial_treatment
table(recur)
yes <- grep("YES",recur)
recur <- rep(0,nrow(LGG.merge))
recur[yes] <- 1
summary(coxph(LGG.os ~ recur)) # Significant (p = 2.28e-07 ***)

# Multivariate survival analysis 
summary(coxph(LGG.os ~ as.numeric(LGG.merge$ENST00000620340.4) + age + recur))# isoform 1  #$coef[1,2] 

summary(coxph(LGG.os ~ as.numeric(LGG.merge$ENST00000262752.4) + age + recur)) # isoform 2 -- not significant 

# Correlation analysis -- whether RSK4 isoform expression is correlated with age and cancer recurrence? 
cor.test(LGG.merge$ENST00000620340.4, LGG.merge$age_at_initial_pathologic_diagnosis) # p-value = 0.0002122, cor = -0.1614438 
ggplot(data = LGG.merge, aes(as.factor(age),ENST00000620340.4, color = as.factor(age))) + geom_boxplot() + 
  stat_compare_means(method = "wilcox.test", paired = FALSE, label.x = 1.25) + # high stage show lower isoform 1 expression
  theme_bw() + 
  labs(x = "Age (>60)", y = "Isoform 1 expression") + 
  scale_color_manual(values = c("blue","red")) +
  theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14),
        plot.title = element_text(size = 16), legend.position = "NONE")

summary(lm(as.numeric(LGG.merge$ENST00000620340.4) ~ as.factor(recur))) # Significant (p= 0.04)
ggplot(data = LGG.merge, aes(as.factor(recur),ENST00000620340.4, color=as.factor(recur))) + geom_boxplot() + 
  stat_compare_means(method = "wilcox.test", paired = FALSE, label.x = 1.25) + # high stage show lower isoform 1 expression
  theme_bw() + 
  labs(x = "Cancer Recurrence", y = "Isoform 1 expression") + 
  scale_color_manual(values = c("blue","red")) +
  theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14),
        plot.title = element_text(size = 16), legend.position = "NONE") 


############################ Adrenocortical cancer (ACC) ############################# 

ACC.rsk <- merge_combine %>% 
  filter(Sample_type == "Adrenocortical Cancer") %>% 
  select(ENST00000620340.4,ENST00000262752.4, ratio, OS, OS.time)

ACC.merge <- merge(ACC.rsk,ACC.pheno,by = 0)

# Create survival object 
ACC.os <- Surv(ACC.merge$OS.time, ACC.merge$OS)

colnames(ACC.merge)
# Covariate 1: Gender 
gender <- as.factor(ACC.merge$gender)
table(gender)
summary(coxph(ACC.os ~ gender)) #Not significant (p = 0.946)

# Covariate 2: Age 
over60 <- which(ACC.merge$age_at_initial_pathologic_diagnosis > 60)
ACC.merge$age_over_60 <- 0 
ACC.merge$age_over_60[over60] <- 1
age <- ACC.merge$age_over_60
summary(coxph(ACC.os ~ age)) # Not significant (p = 0.246)


# Covariate 3: Pathologic stage
stage <- ACC.merge$pathologic_stage
table(stage)
x3<- grep("III",stage)
x4<- grep("IV",stage)
stage.high<-rep(0,nrow(ACC.merge))
stage.high[c(x3,x4)]<-1
summary(coxph(ACC.os ~ stage.high)) # Significant (p = 3.5e-05 ***)

# Covariate 4: cancer recurrence 
recur <- ACC.merge$new_tumor_event_after_initial_treatment
table(recur)
yes <- grep("YES",recur)
recur <- rep(0,nrow(ACC.merge))
recur[yes] <- 1
summary(coxph(ACC.os ~ recur)) # Significant (p = 0.00018 ***)

# Multivariate survival analysis 
summary(coxph(ACC.os ~ as.numeric(ACC.merge$ENST00000620340.4) + stage.high + recur)) # isoform 1 
summary(coxph(ACC.os ~ as.numeric(ACC.merge$ENST00000262752.4) + stage.high + recur)) # isoform 2

# Correlation test -- whether RSK4 isoform expression is correlated with pathologic stage and cancer recurrence in ACC? 
summary(lm(as.numeric(ACC.merge$ENST00000620340.4) ~ as.factor(stage.high)))  # Not significant (p= 0.167)
ggplot(ACC.merge, aes(as.factor(stage.high),ENST00000620340.4 )) + geom_boxplot() # high stage show slightly lower isoform 1 expression

summary(lm(as.numeric(ACC.merge$ENST00000620340.4) ~ as.factor(recur)))  # Not significant (p= 0.194)
ggplot(ACC.merge, aes(as.factor(recur), ENST00000620340.4 )) + geom_boxplot()


############################ Stomach Adenocarcinoma (STAD) ############################# 
STAD.rsk <- merge_combine %>% 
  filter(Sample_type == "Stomach Adenocarcinoma") %>% 
  select(ENST00000620340.4,ENST00000262752.4, ratio, OS, OS.time)

STAD.merge <- merge(STAD.rsk,STAD.pheno,by = 0)

# Create survival object 
STAD.os <- Surv(STAD.merge$OS.time, STAD.merge$OS)

colnames(STAD.merge)
# Covariate 1: Gender 
gender <- as.factor(STAD.merge$gender)
table(gender)
summary(coxph(STAD.os ~ gender)) #Not significant (p = 0.312)

# Covariate 2: Age 
over60 <- which(STAD.merge$age_at_initial_pathologic_diagnosis > 60)
STAD.merge$age_over_60 <- 0 
STAD.merge$age_over_60[over60] <- 1
age <- STAD.merge$age_over_60
summary(coxph(STAD.os ~ age)) # Significant (p =6.67e-05 ***)

# Covariate 3: Pathologic stage
stage <- STAD.merge$pathologic_stage
table(stage)
x3<- grep("III",stage)
x4<- grep("IV",stage)
stage.high<-rep(0,nrow(STAD.merge))
stage.high[c(x3,x4)]<-1
summary(coxph(STAD.os ~ stage.high)) # Significant (p = 0.000234 ***)

# Covariate 4: cancer recurrence 
recur <- STAD.merge$new_tumor_event_after_initial_treatment
table(recur)
yes <- grep("YES",recur)
recur <- rep(0,nrow(STAD.merge))
recur[yes] <- 1
summary(coxph(STAD.os ~ recur)) # Significant (p = 1.35e-08 ***)

# Covariate 5: Metastasis 
meta <- STAD.merge$pathologic_M
table(meta)
M1 <- grep("M1",meta)
metastasis <- rep(0,nrow(STAD.merge))
metastasis[M1] <- 1 
summary(coxph(STAD.os ~ metastasis)) # Significant (p = 0.00198 **)

# Multivariate survival analysis 
summary(coxph(STAD.os ~ as.numeric(STAD.merge$ENST00000620340.4) + age + stage.high + recur + metastasis)) # isoform 1 
summary(coxph(STAD.os ~ as.numeric(STAD.merge$ENST00000262752.4) + age + stage.high + recur + metastasis)) # isoform 2 (not significant)

# Correlation test -- whether RSK4 isoform expression is correlated with pathologic stage and cancer recurrence in ACC? 
cor.test(STAD.merge$ENST00000620340.4, STAD.merge$age_at_initial_pathologic_diagnosis) # p-value = 0.04679, cor = -0.09837304 
ggplot(data = STAD.merge, aes(as.factor(age), ENST00000620340.4)) + geom_boxplot() + stat_compare_means(method = "wilcox.test", paired = FALSE)

summary(lm(as.numeric(STAD.merge$ENST00000620340.4) ~ as.factor(stage.high)))  # Not significant (p= 0.192)
ggplot(STAD.merge, aes(as.factor(stage.high),ENST00000620340.4 )) + geom_boxplot() + stat_compare_means(method = "wilcox.test", paired = FALSE)# high stage show slightly higher isoform 1 expression

summary(lm(as.numeric(STAD.merge$ENST00000620340.4) ~ as.factor(recur)))  # Not significant (p= 0.464)
ggplot(STAD.merge, aes(as.factor(recur), ENST00000620340.4 )) + geom_boxplot() + stat_compare_means(method = "wilcox.test", paired = FALSE)

summary(lm(as.numeric(STAD.merge$ENST00000620340.4) ~ as.factor(metastasis)))  # Not significant (p= 0.98)
ggplot(STAD.merge, aes(as.factor(metastasis), ENST00000620340.4 )) + geom_boxplot() + stat_compare_means(method = "wilcox.test", paired = FALSE)

# isoform 2 
cor.test(STAD.merge$ENST00000262752.4, STAD.merge$age_at_initial_pathologic_diagnosis) # Not significant 
ggplot(data = STAD.merge, aes(age_at_initial_pathologic_diagnosis,ENST00000620340.4)) + geom_point() + geom_smooth(method = "lm")

summary(lm(as.numeric(STAD.merge$ENST00000262752.4) ~ as.factor(stage.high)))  # Not significant (p= 0.336)
ggplot(STAD.merge, aes(as.factor(stage.high),ENST00000262752.4 )) + geom_boxplot() # high stage show slightly higher isoform 1 expression

summary(lm(as.numeric(STAD.merge$ENST00000262752.4) ~ as.factor(recur)))  # Not significant (p= 0.611)
ggplot(STAD.merge, aes(as.factor(recur), ENST00000262752.4 )) + geom_boxplot()

summary(lm(as.numeric(STAD.merge$ENST00000262752.4) ~ as.factor(metastasis)))  # Not significant (p= 0.835)
ggplot(STAD.merge, aes(as.factor(metastasis), ENST00000262752.4 )) + geom_boxplot()


############################ Kidney Clear Cell Carcinoma (KIRC) ############################# 
KIRC.rsk <- merge_combine %>% 
  filter(Sample_type == "Kidney Clear Cell Carcinoma") %>% 
  select(ENST00000620340.4,ENST00000262752.4, ratio, OS, OS.time)

KIRC.merge <- merge(KIRC.rsk,KIRC.pheno,by = 0)

# Create survival object 
KIRC.os <- Surv(KIRC.merge$OS.time, KIRC.merge$OS)

colnames(KIRC.merge)
# Covariate 1: Gender 
gender <- as.factor(KIRC.merge$gender)
table(gender)
summary(coxph(KIRC.os ~ gender)) #Not significant (p = 0.742)

# Covariate 2: Age 
over60 <- which(KIRC.merge$age_at_initial_pathologic_diagnosis > 60)
KIRC.merge$age_over_60 <- 0 
KIRC.merge$age_over_60[over60] <- 1
age <- KIRC.merge$age_over_60
summary(coxph(KIRC.os ~ age)) # Significant (p = 0.000128 ***)

# Covariate 3: Pathologic stage
stage <- KIRC.merge$pathologic_stage
table(stage)
x3<- grep("III",stage)
x4<- grep("IV",stage)
stage.high<-rep(0,nrow(KIRC.merge))
stage.high[c(x3,x4)]<-1
summary(coxph(KIRC.os ~ stage.high)) # Significant (p <2e-16 ***)

# Covariate 4: cancer recurrence 
recur <- KIRC.merge$new_tumor_event_after_initial_treatment
table(recur)
yes <- grep("YES",recur)
recur <- rep(0,nrow(KIRC.merge))
recur[yes] <- 1
summary(coxph(KIRC.os ~ recur)) # Significant (p = 0.017 *)

# Covariate 5: Metastasis 
meta <- KIRC.merge$pathologic_M
table(meta)
M1 <- grep("M1",meta)
metastasis <- rep(0,nrow(KIRC.merge))
metastasis[M1] <- 1 
summary(coxph(KIRC.os ~ metastasis)) # Significant (p <2e-16 ***)

# Multivariate survival analysis 
summary(coxph(KIRC.os ~ as.numeric(KIRC.merge$ENST00000620340.4) + age + stage.high + recur + metastasis)) # isoform 1 
summary(coxph(KIRC.os ~ as.numeric(KIRC.merge$ENST00000262752.4) + age + stage.high + recur + metastasis)) # isoform 2

# Correlation test -- whether RSK4 isoform expression is correlated with pathologic stage and cancer recurrence in ACC? 
cor.test(KIRC.merge$ENST00000620340.4, KIRC.merge$age_at_initial_pathologic_diagnosis) # p-value = 0.02015, cor = -0.1008127 
ggplot(data = KIRC.merge, aes(as.factor(age),ENST00000620340.4, color = as.factor(age))) + geom_boxplot() + 
  stat_compare_means(method = "wilcox.test", paired = FALSE, label.x = 1.25) + 
  theme_bw() + 
  labs(x = "Age (>60)", y = "Isoform 1 expression") + 
  scale_color_manual(values = c("blue","red")) +
  theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14),
                     axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14),
                     plot.title = element_text(size = 16), legend.position = "NONE")


summary(lm(as.numeric(KIRC.merge$ENST00000620340.4) ~ as.factor(stage.high)))  # Significant (p= 0.000381 ***, cor = -0.86)
ggplot(KIRC.merge, aes(as.factor(stage.high),ENST00000620340.4, color = as.factor(stage.high))) + geom_boxplot() + 
  stat_compare_means(method = "wilcox.test", paired = FALSE, label.x = 1.25) + # high stage show lower isoform 1 expression
  theme_bw() + 
  labs(x = "Stage.high", y = "Isoform 1 expression") + 
  scale_color_manual(values = c("blue","red")) +
  theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14),
        plot.title = element_text(size = 16), legend.position = "NONE")


summary(lm(as.numeric(KIRC.merge$ENST00000620340.4) ~ as.factor(recur)))  # Significant (p= 0.0412 *, cor = -1.1652)
ggplot(KIRC.merge, aes(as.factor(recur), ENST00000620340.4 )) + geom_boxplot()+ stat_compare_means(method = "wilcox.test", paired = FALSE)

summary(lm(as.numeric(KIRC.merge$ENST00000620340.4) ~ as.factor(metastasis)))  # Significant (p= 0.00343 **, cor = -0.9832)
ggplot(KIRC.merge, aes(as.factor(metastasis), ENST00000620340.4 )) + geom_boxplot() + stat_compare_means(method = "wilcox.test", paired = FALSE)

#isoform 2
summary(lm(as.numeric(KIRC.merge$ENST00000262752.4) ~ as.factor(stage.high)))  # Significant (p= 0.0024 **  , cor = -1.0795)
ggplot(KIRC.merge, aes(as.factor(stage.high),ENST00000262752.4, color = as.factor(stage.high) )) + geom_boxplot()+ 
  stat_compare_means(method = "wilcox.test", paired = FALSE, label.x = 1.25) + # high stage show lower isoform 2 expression
  theme_bw() + 
  labs(x = "Stage.high", y = "Isoform 2 expression") + 
  scale_color_manual(values = c("blue","red")) +
  theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14),
        plot.title = element_text(size = 16), legend.position = "NONE")

summary(lm(as.numeric(KIRC.merge$ENST00000262752.4) ~ as.factor(metastasis)))
ggplot(KIRC.merge, aes(as.factor(metastasis),ENST00000262752.4, color = as.factor(metastasis))) + geom_boxplot()+ 
  stat_compare_means(method = "wilcox.test", paired = FALSE, label.x = 1.25) + # high stage show lower isoform 2 expression
  theme_bw() + 
  labs(x = "Metastasis", y = "Isoform 2 expression") + 
  scale_color_manual(values = c("blue","red")) +
  theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14),
        plot.title = element_text(size = 16), legend.position = "NONE")

############################ Cervical & Endocervical Cancer (CESC) ############################# 
CESC.rsk <- merge_combine %>% 
  filter(Sample_type == "Cervical & Endocervical Cancer" ) %>% 
  select(ENST00000620340.4,ENST00000262752.4, ratio, OS, OS.time)

CESC.merge <- merge(CESC.rsk,CESC.pheno,by = 0)

# Create survival object 
CESC.os <- Surv(CESC.merge$OS.time, CESC.merge$OS)

colnames(CESC.merge)
# Covariate 1: Gender 
gender <- as.factor(CESC.merge$gender)
table(gender)  ## all patients are females 

# Covariate 2: Age 
over60 <- which(CESC.merge$age_at_initial_pathologic_diagnosis > 60)
CESC.merge$age_over_60 <- 0 
CESC.merge$age_over_60[over60] <- 1
age <- CESC.merge$age_over_60
summary(coxph(CESC.os ~ age)) # Significant (p = 0.0196 *)

# Covariate 3: Clinical stage
stage <- CESC.merge$clinical_stage
table(stage)
x3<- grep("III",stage)
x4<- grep("IV",stage)
stage.high<-rep(0,nrow(CESC.merge))
stage.high[c(x3,x4)]<-1
summary(coxph(CESC.os ~ stage.high)) # Significant (p = 0.000371 ***)

# Covariate 4: Cancer recurrence 
recur <- CESC.merge$new_tumor_event_after_initial_treatment
table(recur)
yes <- grep("YES",recur)
recur <- rep(0,nrow(CESC.merge))
recur[yes] <- 1
summary(coxph(CESC.os ~ recur)) # Significant (p = 6.77e-11 ***)

# Covariate 5: Metastasis 
meta <- CESC.merge$pathologic_M
table(meta)
M1 <- grep("M1",meta)
metastasis <- rep(0,nrow(CESC.merge))
metastasis[M1] <- 1 
summary(coxph(CESC.os ~ metastasis)) # Not significant (p = 0.11)

# Multivariate survival analysis 
summary(coxph(CESC.os ~ as.numeric(CESC.merge$ENST00000262752.4) + age + stage.high + recur )) # isoform 2

# Correlation test -- whether RSK4 isoform expression is correlated with pathologic stage and cancer recurrence in ACC? 
summary(lm(as.numeric(CESC.merge$ENST00000262752.4) ~ as.factor(stage.high)))  # Not significant (p= 0.345)
ggplot(CESC.merge, aes(as.factor(stage.high),ENST00000262752.4 )) + geom_boxplot() + stat_compare_means(method = "wilcox.test", paired = FALSE)# high stage show slightly higher isoform 1 expression

summary(lm(as.numeric(CESC.merge$ENST00000262752.4) ~ as.factor(recur)))  # Not significant (p= 0.722) isoform 1 # No significant (p= 0.593) isoform 2 
ggplot(CESC.merge, aes(as.factor(recur), ENST00000262752.4 )) + geom_boxplot() + stat_compare_means(method = "wilcox.test", paired = FALSE)

summary(lm(as.numeric(CESC.merge$ratio) ~ as.factor(stage.high)))  # Not significant (p= 0.059)
ggplot(CESC.merge, aes(as.factor(stage.high),ratio )) + geom_boxplot() + stat_compare_means(method = "wilcox.test", paired = FALSE)# high stage show slightly higher isoform 1 expression

############################ Rectum Adenocarcinoma (READ) ############################# 
READ.rsk <- merge_combine %>% 
  filter(Sample_type == "Rectum Adenocarcinoma" ) %>% 
  select(ENST00000620340.4,ENST00000262752.4, ratio, OS, OS.time)

READ.merge <- merge(READ.rsk,READ.pheno,by = 0)

# Create survival object 
READ.os <- Surv(READ.merge$OS.time, READ.merge$OS)

colnames(READ.merge)
# Covariate 1: Gender 
gender <- as.factor(READ.merge$gender)
table(gender)  
summary(coxph(READ.os ~ gender)) # Not significant (p=0.361)

# Covariate 2: Age 
age <- READ.merge$age_at_initial_pathologic_diagnosis
table(age)
summary(coxph(READ.os ~ age)) # Significant (p = 0.000707 ***)

# Covariate 3: Clinical stage
stage <- READ.merge$pathologic_stage
table(stage)
x3<- grep("III",stage)
x4<- grep("IV",stage)
stage.high<-rep(0,nrow(READ.merge))
stage.high[c(x3,x4)]<-1
summary(coxph(READ.os ~ stage.high)) # Not significant (p = 0.269)

# Covariate 4: Cancer recurrence 
recur <- READ.merge$new_tumor_event_after_initial_treatment
table(recur)
yes <- grep("YES",recur)
recur <- rep(0,nrow(READ.merge))
recur[yes] <- 1
summary(coxph(READ.os ~ recur)) # Not significant (p = 0.104)

# Covariate 5: Metastasis 
meta <- READ.merge$pathologic_M
table(meta)
M1 <- grep("M1",meta)
metastasis <- rep(0,nrow(READ.merge))
metastasis[M1] <- 1 
summary(coxph(READ.os ~ metastasis)) # Not significant (p = 0.766)


# Multivariate survival analysis 
summary(coxph(READ.os ~ as.numeric(READ.merge$ENST00000262752.4) + age)) # isoform 1 

# Correlation test -- whether RSK4 isoform expression is correlated with pathologic stage and cancer recurrence in ACC? 
cor.test(READ.merge$ENST00000262752.4, age) # Not significant (p= 0.7244)
ggplot(READ.merge, aes(age_at_initial_pathologic_diagnosis, ENST00000262752.4)) + geom_point() + geom_smooth(method = "lm")# high stage show slightly higher isoform 1 expression




