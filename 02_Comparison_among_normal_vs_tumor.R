##　Project Aim: Differential expression analysis of RSK4 isoforms in tumour and normal tissue samples using TCGA and GTEx data

# Objectives:
	# 1. Examine the data normality using density plot, Q-Q plot and Shapiro test before comparison
	# 2. Compare the RSK4 isoform expression among normal tissues using GTEx samples 
	# 3. Compare the RSK4 isoform expression among tumour samples using TCGA samples 
	# 4. Compare the RSK4 isoform expression between normal tissues and corresponding cancerous tissues (GTEx vs TCGA)

# Load packages 
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(FSA)  # Dunn test 
library(GmAMisc) # KwPlot 
library(rstatix) # Cohen test
library(effsize)
library(EnvStats) # Outlier detection - Rosner test
library(viridis) 

# Load data 
setwd("/Users/chensisi/Documents/RNAseq/")
load("RSK4_data.RData")


## Objective 1: Data normality test ##

	# 1. Density plot #

	# TCGA 
ggplot(TCGA_by_sample, aes(x=ENST00000620340.4, fill = Sample_type)) + 
  geom_density(alpha=.3) + 
  labs(title = "Density plot for RSK4 isoform 1 expression in TCGA samples")

ggplot(TCGA_by_sample, aes(x=ENST00000262752.4, fill = Sample_type)) + 
  geom_density(alpha=.3) + 
  labs(title = "Density plot for RSK4 isoform 2 expression in TCGA samples")


	# GTEX
ggplot(GTEX_by_sample, aes(x=ENST00000620340.4, fill = Sample_type)) + 
  geom_density(alpha=.3) + 
  labs(title = "Density plot for RSK4 isoform 1 expression in GTEX samples")

ggplot(GTEX_by_sample, aes(x=ENST00000262752.4, fill = Sample_type)) + 
  geom_density(alpha=.3) + 
  labs(title = "Density plot for RSK4 isoform 2 expression in GTEX samples")

	# 2. Q-Q plot #

  # By project: TCGA vs GTEX 
qplot(sample = ENST00000620340.4, data = merge1 , color= Project) + 
  labs(title = "Q-Q plot for RSK4 isoform 1 expression", 
       y = "log2(x+0.0001) transfromed TPM", x = "Threoretical") 

qplot(sample = ENST00000262752.4, data = merge1 , color= Project) + 
  labs(title = "Q-Q plot for RSK4 isoform 1 expression", 
       y = "log2(x+0.0001) transfromed TPM", x = "Threoretical")
  # TCGA 
qplot(sample = ENST00000620340.4, data = TCGA_by_sample, color = Sample_type) + 
  labs(title = "Q-Q plot for RSK4 isoform 1 expression in TCGA samples", 
       y = "log2(x+0.0001) transformed TPM", x = "Threoretical")

qplot(sample = ENST00000262752.4, data = TCGA_by_sample, color = Sample_type) + 
  labs(title = "Q-Q plot for RSK4 isoform 2 expression in TCGA samples", 
     y = "log2(x+0.0001) transformed TPM", x = "Threoretical")

  # GTEX
qplot(sample = ENST00000620340.4, data = GTEX_by_sample, color = Sample_type) + 
  labs(title = "Q-Q plot for RSK4 isoform 1 expression in GTEX samples", 
       y = "log2(x+0.0001) transformed TPM", x = "Threoretical")

qplot(sample = ENST00000262752.4, data = GTEX_by_sample, color = Sample_type) + 
  labs(title = "Q-Q plot for RSK4 isoform 2 expression in GTEX samples", 
       y = "log2(x+0.0001) transformed TPM", x = "Threoretical")

	# 3. Shapiro test #
 		#TCGA by sample types
tcga_normality <- data.frame()

TCGA_samples <- unique(TCGA_by_sample$Sample_type)

for (x in TCGA_samples){
  a <- TCGA_by_sample %>%
    select(ENST00000620340.4,ENST00000262752.4,Sample_type) %>%
    filter(Sample_type == x)
  shapiro_p_value <- shapiro.test(a$ENST00000620340.4)[[2]] #ENST00000262752.4
  result <- data.frame(x, shapiro_p_value)
  tcga_normality <- rbind(tcga_normality, result)
}

tcga_normality %>% filter(shapiro_p_value >= 0.05) # Only results with p-value larger than 0.05 are considered as normally distributed 


		# GTEX by sample types
gtex_normality <- data.frame()

GTEX_samples <- unique(GTEX_by_sample$Sample_type)

for (x in GTEX_samples){
  a <- GTEX_by_sample %>%
    select(ENST00000620340.4,ENST00000262752.4,Sample_type) %>% # ENST00000262752.4
    filter(Sample_type == x)
  shapiro_p_value <- shapiro.test(a$ENST00000620340.4)[[2]]
  result <- data.frame(x, shapiro_p_value)
  gtex_normality <- rbind(gtex_normality, result) #  Only results with p-value larger than 0.05 are considered as normally distributed 
}


## Objective 2: RSK4 isoform expression comparison among normal tissues ##
	# Kruskal Wallis test 
iso1.kw <- kruskal.test(ENST00000620340 ~ Sample_type, data = GTEX_by_sample) # Isoform 1 expression 
iso2.kw <- kruskal.test(ENST00000262752 ~ Sample_type, data = GTEX_by_sample) # Isoform 2 expression 

	# Multiple pairwise comparison -- Dunn test 
iso1.d <- dunnTest(ENST00000620340.4 ~ as.factor(Sample_type), data = GTEX_by_sample, method="bonferroni")$res
iso2.d <- dunnTest(ENST00000262752.4 ~ as.factor(Sample_type), data = GTEX_by_sample, method="bonferroni")$res

iso1.d %>%  
	filter(P.adj < 0.05 ) %>%  # retain only significant results 
	arrange(P.adj)

iso2.d %>% 
	filter(P.adj < 0.05 ) %>%
	arrange(P.adj)
	

## Objective 3: Compare the RSK4 isoform expression among tumour samples using TCGA samples 
	# Kruskal Wallis test 
iso1.kw <- kruskal.test(ENST00000620340 ~ Sample_type, data = TCGA_by_sample) # Isoform 1 expression 
iso2.kw <- kruskal.test(ENST00000262752 ~ Sample_type, data = TCGA_by_sample) # Isoform 2 expression 

	# Multiple pairwise comparison -- Dunn test 
iso1.d <- dunnTest(ENST00000620340.4 ~ as.factor(Sample_type), data = TCGA_by_sample, method="bonferroni")$res
iso2.d <- dunnTest(ENST00000262752.4 ~ as.factor(Sample_type), data = TCGA_by_sample, method="bonferroni")$res

iso1.d %>%  
	filter(P.adj < 0.05 ) %>%  # retain only significant results 
	arrange(P.adj)

iso2.d %>% 
	filter(P.adj < 0.05 ) %>%
	arrange(P.adj)


## Objective 4: Compare the RSK4 isoform expression between normal tissues and corresponding cancerous tissues (GTEx vs TCGA)

# Compare the RSK4 isoform expression between TCGA & GTEX samples in general 
wilcox_results1 <- wilcox.test(ENST00000620340.4 ~ Project, data = merge1)
wilcox_results1 

wilcox_results2 <- wilcox.test(ENST00000262752.4 ~ Project, data = merge1)
wilcox_results2

# Data visualization -- compare the isoform 1 and 2 expression in general between GTEX and TCGA data 
    # Isoform 1 
merge1 %>% 
  ggplot(aes(Project,ENST00000620340.4, color = Project)) + geom_boxplot() + 
  stat_compare_means(method ="wilcox.test", label.y = 4, label.x = 1.25, size =4.5) + 
  my_theme + blank_theme + 
  scale_color_manual(labels = c("GTEX","TCGA"), values =c("#00BFC4","#F8766D"))

  # Isoform 2
merge1 %>% 
  ggplot(aes(Project,ENST00000262752.4, color = Project)) + geom_boxplot() + 
  stat_compare_means(method ="wilcox.test", label.y = 3.5, label.x = 1.25, size =4.5) + 
  my_theme + blank_theme + 
  scale_color_manual(labels = c("GTEX","TCGA"), values =c("#00BFC4","#F8766D"))

merge1 %>% 
  ggplot(aes(Project,ENST00000262752.4, color = Project)) + geom_boxplot() + 
  theme_light() + my_theme + theme(legend.position = "top", legend.background = element_rect(fill = "white")) + 
  scale_color_manual(labels = c("GTEX","TCGA"), values =c("#00BFC4","#F8766D"))


# 1) Create dataframe (normal tissue with corresponding cancer types)
# Bladder
bladder<- merge1 %>% 
  group_by(Project) %>% 
  select(Project, ENST00000620340.4, ENST00000262752.4,  Sample_type) %>%
  filter(Project == "TCGA" | Project == "GTEX") %>% 
  filter(Sample_type == "Bladder" | Sample_type == "Bladder Urothelial Carcinoma") 

# Breast 
breast <- merge1 %>% 
  group_by(Project) %>% 
  select(Project, ENST00000620340.4, ENST00000262752.4,  Sample_type) %>%
  filter(Project == "TCGA" | Project == "GTEX") %>% 
  filter(Sample_type == "Breast" | Sample_type == "Breast Invasive Carcinoma") 

# Cervix & uterus 
cerv_uter <- merge1 %>% 
  group_by(Project) %>% 
  select(Project, ENST00000620340.4, ENST00000262752.4, Sample_type) %>%
  filter(Project == "TCGA" | Project == "GTEX") %>% 
  filter(Sample_type == "Cervix Uteri" | Sample_type == "Uterus"|
           Sample_type  == "Cervical & Endocervical Cancer" | Sample_type == "Uterine Carcinosarcoma") #Sample_type == "Uterine Corpus Endometrioid Carcinoma"

# Brain 
brain <- merge1 %>% 
  group_by(Project) %>% 
  select(Project, ENST00000620340.4, ENST00000262752.4, Sample_type) %>%
  filter(Project == "TCGA" | Project == "GTEX") %>% 
  filter(Sample_type == "Brain" | Sample_type == "Glioblastoma Multiforme"|
           Sample_type  == "Brain Lower Grade Glioma" ) 

# Fallopian Tube 
fallopian <- merge1 %>% 
  group_by(Project) %>% 
  select(Project, ENST00000620340.4, ENST00000262752.4, Sample_type) %>%
  filter(Project == "TCGA" | Project == "GTEX") %>% 
  filter(Sample_type == "Fallopian Tube" | Sample_type == "Ovarian Serous Cystadenocarcinoma") 

# Lung 
lung <- merge1 %>% 
  group_by(Project) %>% 
  select(Project, ENST00000620340.4, ENST00000262752.4, Sample_type) %>%
  filter(Project == "TCGA" | Project == "GTEX") %>% 
  filter(Sample_type == "Lung" | Sample_type == "Lung Adenocarcinoma" | 
           Sample_type == "Lung Squamous Cell Carcinoma") 

# Prostate 
prostate <- merge1 %>% 
  group_by(Project) %>% 
  select(Project, ENST00000620340.4, ENST00000262752.4, Sample_type) %>%
  filter(Project == "TCGA" | Project == "GTEX") %>% 
  filter(Sample_type == "Prostate" | Sample_type == "Prostate Adenocarcinoma") 

# Testis
testis <- merge1 %>% 
  group_by(Project) %>% 
  select(Project, ENST00000620340.4, ENST00000262752.4, Sample_type) %>%
  filter(Project == "TCGA" | Project == "GTEX") %>% 
  filter(Sample_type == "Testis" | Sample_type == "Testicular Germ Cell Tumor") 

# Esophagus
esophagus <- merge1 %>% 
  group_by(Project) %>% 
  select(Project, ENST00000620340.4, ENST00000262752.4, Sample_type) %>%
  filter(Project == "TCGA" | Project == "GTEX") %>% 
  filter(Sample_type == "Esophagus" | Sample_type == "Esophageal Carcinoma") 

# Pancreas
pancreas <- merge1 %>% 
  group_by(Project) %>% 
  select(Project, ENST00000620340.4, ENST00000262752.4, Sample_type) %>%
  filter(Project == "TCGA" | Project == "GTEX") %>% 
  filter(Sample_type == "Pancreas" | Sample_type == "Pancreatic Adenocarcinoma") 

# Kidney
kidney <- merge1 %>% 
  group_by(Project) %>% 
  select(Project, ENST00000620340.4, ENST00000262752.4, Sample_type) %>%
  filter(Project == "TCGA" | Project == "GTEX") %>% 
  filter(Sample_type == "Kidney" | Sample_type == "Kidney Papillary Cell Carcinoma"|
           Sample_type == "Kidney Clear Cell Carcinoma"| Sample_type == "Kidney Chromophobe") 

# Liver
liver <- merge1 %>% 
  group_by(Project) %>% 
  select(Project, ENST00000620340.4, ENST00000262752.4, Sample_type) %>%
  filter(Project == "TCGA" | Project == "GTEX") %>% 
  filter(Sample_type == "Liver" | Sample_type == "Liver Hepatocellular Carcinoma") 

# Muscle 
muscle <- merge1 %>% 
  group_by(Project) %>% 
  select(Project, ENST00000620340.4, ENST00000262752.4, Sample_type) %>%
  filter(Project == "TCGA" | Project == "GTEX") %>% 
  filter(Sample_type == "Muscle"| Sample_type == "Sarcoma" )  

# Colon 
colon <- merge1 %>% 
  group_by(Project) %>% 
  select(Project, ENST00000620340.4, ENST00000262752.4, Sample_type) %>%
  filter(Project == "TCGA" | Project == "GTEX") %>% 
  filter(Sample_type == "Colon" | Sample_type == "Colon Adenocarcinoma")  

# Stomach 
stomach <- merge1 %>% 
  group_by(Project) %>% 
  select(Project, ENST00000620340.4, ENST00000262752.4, Sample_type) %>%
  filter(Project == "TCGA" | Project == "GTEX") %>% 
  filter(Sample_type == "Stomach" | Sample_type == "Stomach Adenocarcinoma") 

# Skin 
skin <- merge1 %>% 
  group_by(Project) %>% 
  select(Project, ENST00000620340.4, ENST00000262752.4, Sample_type) %>%
  filter(Project == "TCGA" | Project == "GTEX") %>% 
  filter(Sample_type == "Skin" | Sample_type == "Skin Cutaneous Melanoma") 

# Thyroid  
thyroid <- merge1 %>% 
  group_by(Project) %>% 
  select(Project, ENST00000620340.4, ENST00000262752.4, Sample_type) %>%
  filter(Project == "TCGA" | Project == "GTEX") %>% 
  filter(Sample_type == "Thyroid" | Sample_type == "Thyroid Carcinoma") 

# Blood   
blood <- merge1 %>% 
  group_by(Project) %>% 
  select(Project, ENST00000620340.4, ENST00000262752.4, Sample_type) %>%
  filter(Project == "TCGA" | Project == "GTEX") %>% 
  filter(Sample_type == "Blood" | Sample_type == "Acute Myeloid Leukemia"|
           Sample_type == "Diffuse Large B-Cell Lymphoma") 

# Adrenal Gland
adrenal <- merge1 %>% 
  group_by(Project) %>% 
  select(Project, ENST00000620340.4, ENST00000262752.4, Sample_type) %>%
  filter(Project == "TCGA" | Project == "GTEX") %>% 
  filter(Sample_type == "Adrenal Gland" | Sample_type == "Adrenocortical Cancer"|
           Sample_type == "Pheochromocytoma & Paraganglioma") 



### Wilcoxon test to compare the expression difference between normal tissues and corresponding cancer types 
# Data visualization 
list.dfs <- list(bladder, breast, cerv_uter, brain, fallopian,
                 lung, prostate, testis, esophagus, pancreas,
                 kidney, liver, muscle, colon, stomach, skin,
                 thyroid, blood, adrenal)

groups <- lapply(list.dfs, function(x) {
                    x$Group <- unique(x$Sample_type[1]) 
                    x} )
tumor_vs_normal <- data.frame()
for (i in 1:length(groups)){
  data <- as.data.frame(groups[[i]])
  tumor_vs_normal <- rbind(tumor_vs_normal, data)
}

# Multi-panel plot -- compare RSK4 isoform expression between TCGA and GTEX samples 
p <- ggplot(tumor_vs_normal, aes(x=Project, y = ENST00000262752.4, color = Project)) +  # ENST00000262752.4/ ENST00000620340.4
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1)) + theme(plot.title = element_text(hjust = 0.5)) + facet_wrap(~ Group)


# Add p-value 
p + stat_compare_means(method ="wilcox.test", label.y = 1.5, label.x = 1.5, label = "p.signif") + 
  theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12),
        legend.text = element_text(size = 10), legend.title = element_text(size =12),
        strip.text = element_text(size = 14, face = "bold")) + theme_linedraw() + 
  scale_color_manual(labels = c("GTEX","TCGA"), values =c("#00BFC4","#F8766D"))


# Statistic tests (Dunn's test + Cohen's test)
	## Dunn's test 
	#1. Compare by actual expression (used normalized TPM value)
iso1.results <- data.frame()
iso2.results <- data.frame()
for (i in 1:length(list.dfs)){
  df <- as.data.frame(list.dfs[[i]])
  isoform1 <- dunn_test(ENST00000620340.4 ~ Sample_type, data = df)
  isoform2 <- dunn_test(ENST00000262752.4 ~ Sample_type, data = df)
  iso1.results <- rbind(iso1.results,isoform1)
  iso2.results <- rbind(iso2.results,isoform2)
}

	#2. Compare by ratio (isoform2/isoform1)
ratio.results <- data.frame()
for (i in 1:length(list.dfs_ratio)){
  df <- as.data.frame(list.dfs_ratio[[i]])
  ratio <- dunn_test(ratio ~ Sample_type, data = df)
  ratio.results <- rbind(ratio.results,ratio)
}

	## Cohen's test 
	#1. Compare by actual expression (used normalized TPM value)
cohen.test_tpm <- function(x) {
  x <- as.data.frame(x)
  x$Sample_type <- as.factor(x$Sample_type)
  results <- x %>% cohens_d(ENST00000620340.4 ~ Sample_type, var.equal = FALSE) # ENST00000262752.4
  return(results)
} # create the function for cohen test, to do cohen test, we need data frame and formular with column as factor 


results <- map(list.dfs, ~cohen.test_tpm(.)) # apply the function on multiple dataframes 

cohen.results_tpm <- data.frame()
for(i in 1:length(results)){
  cohen.results_tpm <- rbind(cohen.results_tpm, results[[i]])
}

colnames(cohen.results_tpm)<-c("TranscriptIDs","Group1","Group2","effsize","no.GTEX","no.TCGA","magnitude")


	#2. Compare by ratio (isoform2/isoform1)
cohen.test_ratio <- function(x) {
  x <- as.data.frame(x)
  x$Sample_type <- as.factor(x$Sample_type)
  results <- x %>% cohens_d(ratio ~ Sample_type, var.equal = FALSE) # IF2
  return(results)
}  

results <- map(list.dfs_ratio, ~cohen.test_ratio(.)) # apply the function on multiple dataframes 

cohen.results.ratio <- data.frame()  
for(i in 1:length(results)){
  cohen.results.ratio <- rbind(cohen.results.ratio, results[[i]])
}

colnames(cohen.results.ratio)<-c("TranscriptIDs","Group1","Group2","effsize","no.GTEX","no.TCGA","magnitude")


### Combine cohen results with Dunn test results and export(by TPM and ratio) ####
	# TPM 
combine.tpm <- cbind(iso1.results, cohen.results_tpm) # iso2.results
combine.tpm <- combine.tpm %>%
  select(!c("TranscriptIDs","group1","group2","no.GTEX","no.TCGA"))%>%
  arrange(desc(abs(effsize)))
combine.tpm[1:5,]

write.table(combine.tpm, "TCGAvsGTEX.iso1.dunn.cohen.merge.csv", sep=",", row.names = FALSE)


	# ratio 
combine.ratio <- cbind(ratio.results, cohen.results.ratio) # IF1.results, cohen.results.IF1
combine.ratio <- combine.ratio %>%
  select(!c("TranscriptIDs","Group1","Group2","no.GTEX","no.TCGA"))%>%
  arrange(desc(abs(effsize)))
combine.ratio[1:5,]

write.table(combine1, "TCGAvsGTEX.ratio.dunn.cohen.merge.csv", sep=",", row.names = FALSE)

