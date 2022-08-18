# Load packages 
library(gridExtra)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(FSA)  # Dunn test 
library(GmAMisc) # KwPlot 
library(rstatix) # Cohen test
library(effsize)
library(survival)
library(survminer)
library(dplyr)
library(tidyverse)

# Import data 
setwd("/Users/chensisi/Documents/RNAseq/")
load("TCGA_GTEX_merge.RData")
load("RSK4_data.RData")

# set the theme for ggplot 
my_theme <- theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14),
                  axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14),
                  plot.title = element_text(size = 16), legend.text = element_text(size = 12), legend.title = element_text(size =12)) + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"),
                  panel.border = element_blank())

blank_theme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text.x = element_blank(), axis.line = element_line(color = "black"),
                     panel.border = element_blank())


## Proportion of 4 transcript isoforms of RSK4 in TCGA and GTEX samples
merge1_denorm <- merge1 %>% 
  mutate(ENST00000620340.4 = 2^ENST00000620340.4,
         ENST00000495332.1 = 2^ENST00000495332.1,
         ENST00000460730.1 = 2^ENST00000460730.1,
         ENST00000262752.4 = 2^ENST00000262752.4)

TCGA_mean <- merge1_denorm %>%
  filter(Project == "TCGA") %>%
  select("ENST00000620340.4","ENST00000495332.1","ENST00000460730.1","ENST00000262752.4","Project","Sample_type") %>%
  group_by(Sample_type) %>%
  summarize(iso1= mean(ENST00000620340.4), iso2 = mean(ENST00000262752.4),
            iso3= mean(ENST00000495332.1), iso4 = mean(ENST00000460730.1)) %>%
  mutate(total= iso1 + iso2 + iso3 + iso4) 

## Get abbreviations for cancer types 
unique(merge3$cancer.type.abbreviation)
abbreviations <- merge3 %>% select(Sample_type, cancer.type.abbreviation) %>% distinct()
rownames(abbreviations) <- NULL
TCGA_mean <- merge(TCGA_mean, abbreviations)

GTEX_mean <- merge1_denorm %>%
  filter(Project == "GTEX") %>%
  select("ENST00000620340.4","ENST00000495332.1","ENST00000460730.1","ENST00000262752.4","Project","Sample_type") %>%
  group_by(Sample_type) %>%
  summarize(iso1= mean(ENST00000620340.4), iso2 = mean(ENST00000262752.4),
            iso3= mean(ENST00000495332.1), iso4 = mean(ENST00000460730.1)) %>%
  mutate(total= iso1 + iso2 + iso3 + iso4)


##############################################################################################################
################################################## Figure 1 ##################################################
######### Fig.1A #######
# Stacked barplot showing the proportion of RSK4 isoforms in each cancer type
p1 <- ggplot(TCGA_mean, aes(x= Sample_type, y= total)) + geom_bar(stat = "identity", width = 0.9) + 
  theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
         axis.ticks.x=element_blank()) + 
  ylab("Total expression")

p2 <- TCGA_mean %>%
  gather(isoform, expression, iso1:iso4) %>%
  ggplot(aes(x = cancer.type.abbreviation, y = expression, fill = isoform)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_y_continuous(labels = scales::percent_format()) + 
  labs(y = "Relative proportion of RSK4 isoform expression", x = "Cancer types") + 
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1), 
        plot.title = element_text(hjust = 0.5), 
        legend.position = "none")  # remove legends 

plot_grid(p1, p2, labels=c("A", "B"), ncol = 1 , align = "v", axis = "bt",rel_widths = c(1,1.3), rel_heights = c(1.4,4))

# Stacked barplot showing the proportion of RSK4 isoforms in each normal tissue 
p1 <- ggplot(GTEX_mean, aes(x= Sample_type, y= total)) + geom_bar(stat = "identity", width = 0.9) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  ylab("Total expression")

p2 <- GTEX_mean %>%
  gather(isoform, expression, iso1:iso4) %>%
  ggplot(aes(x = Sample_type, y = expression, fill = isoform)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_y_continuous(labels = scales::percent_format()) + 
  labs(y = "Relative proportion of RSK4 isoform expression", x = "Normal Tissues") + 
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1), 
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")  # remove legends 

plot_grid(p1, p2, labels=c("A", "B"), ncol = 1 , align = "v", axis = "bt",rel_widths = c(1,1.3), rel_heights = c(1.4,4))

######### Fig.1B #######
## Boxplot showing the distribution of RSK4 isoform expression among normal and tumor samples
  # Categorize GTEX and TCGA samples into high, medium and low isoform expression groups ####
    #TCGA 
TCGA_by_sample_median <- TCGA_by_sample %>% 
  group_by(Sample_type) %>%
  mutate(median1 = median(ENST00000620340.4), median2 = median(ENST00000262752.4))

TCGA_by_sample_median$Expression <- NA

iso1lower_bound <- quantile(TCGA_by_sample$ENST00000620340.4, 0.25)
iso2lower_bound <- quantile(TCGA_by_sample$ENST00000262752.4, 0.25)
iso1upper_bound <- quantile(TCGA_by_sample$ENST00000620340.4, 0.75)
iso2upper_bound <- quantile(TCGA_by_sample$ENST00000262752.4, 0.75)

#High 
iso1.high <- which(TCGA_by_sample_median$median1 > iso1upper_bound)
iso2.high <- which(TCGA_by_sample_median$median2 > iso2upper_bound)

#Medium 
iso1.medium <- which(TCGA_by_sample_median$median1 > iso1lower_bound & TCGA_by_sample_median$median1 < iso1upper_bound)
iso2.medium <- which(TCGA_by_sample_median$median2 > iso2lower_bound & TCGA_by_sample_median$median2 < iso2upper_bound)

#Low 
iso1.low <- which(TCGA_by_sample_median$median1 < iso1lower_bound | TCGA_by_sample_median$median1 == iso1lower_bound)
iso2.low <- which(TCGA_by_sample_median$median2 < iso2lower_bound | TCGA_by_sample_median$median2 == iso2lower_bound)

#Plot the figure for TCGA 
## Isoform 1 
TCGA_by_sample_median$Expression[iso1.high] <- "High"
TCGA_by_sample_median$Expression[iso1.medium] <- "Medium"
TCGA_by_sample_median$Expression[iso1.low] <- "Low"
TCGA_by_sample_median <- TCGA_by_sample_median %>%
  arrange(desc(median1))

ggplot(TCGA_by_sample_median, aes(x= ENST00000620340.4, y = reorder(Sample_type,median1), color= Expression)) + 
  geom_boxplot() + # cancer.type.abbreviation
  labs( x = "Expression (log2 normalized TPM)", y = "TCGA",title = "Isoform 1") + 
  theme(axis.text.y = element_text(size = 12), axis.title = element_text(size = 14), 
        plot.title = element_text( face = "bold", hjust = 0.5)) + 
  geom_vline(xintercept= c(iso1upper_bound, iso1lower_bound) , linetype="dashed", color = "red")

## Isoform 2
TCGA_by_sample_median$Expression <- NA 
TCGA_by_sample_median$Expression[iso2.high] <- "High"
TCGA_by_sample_median$Expression[iso2.medium] <- "Medium"
TCGA_by_sample_median$Expression[iso2.low] <- "Low"
TCGA_by_sample_median <- TCGA_by_sample_median %>%
  arrange(desc(median2))

ggplot(TCGA_by_sample_median, aes(x= ENST00000262752.4, y = reorder(Sample_type,median2), color= Expression)) + 
  geom_boxplot() + #cancer.type.abbreviation
  labs( x = "Expression (log2 normalized TPM)", y = "TCGA",title = "Isoform 2") + 
  theme(axis.text.y = element_text(size = 12), axis.title = element_text(size = 14), 
        plot.title = element_text( face = "bold", hjust = 0.5)) + 
  geom_vline(xintercept= c(iso1upper_bound, iso1lower_bound) , linetype="dashed", color = "red")

    #GTEX 
GTEX_by_sample_median <- GTEX_by_sample %>% 
  group_by(Sample_type) %>%
  mutate(median1 = median(ENST00000620340.4), median2 = median(ENST00000262752.4)) 

GTEX_by_sample_median$Expression <- NA

iso1lower_bound <- quantile(GTEX_by_sample$ENST00000620340.4, 0.25)
iso2lower_bound <- quantile(GTEX_by_sample$ENST00000262752.4, 0.25)
iso1upper_bound <- quantile(GTEX_by_sample$ENST00000620340.4, 0.75)
iso2upper_bound <- quantile(GTEX_by_sample$ENST00000262752.4, 0.75)

# Isoform 1 
high <- which(GTEX_by_sample_median$median1 >= iso1upper_bound)
medium <- which(GTEX_by_sample_median$median1 < iso1upper_bound & GTEX_by_sample_median$median1 >= iso1lower_bound)
low <- which(GTEX_by_sample_median$median1 < iso1lower_bound )

# Isoform 2 
high <- which(GTEX_by_sample_median$median2 >= iso2upper_bound)
medium <- which(GTEX_by_sample_median$median2 < iso2upper_bound & GTEX_by_sample_median$median2 >= iso2lower_bound)
low <- which(GTEX_by_sample_median$median2 < iso2lower_bound)
#low <- which(GTEX_by_sample_median$median2 < -9.96)

GTEX_by_sample_median$Expression[high] <- "High"
GTEX_by_sample_median$Expression[medium] <- "Medium"
GTEX_by_sample_median$Expression[low] <- "Low"
GTEX_by_sample_median <- GTEX_by_sample_median %>%
  arrange(desc(median2))

# Plot the figure  
ggplot(GTEX_by_sample_median, aes(x= ENST00000262752.4, y = reorder(Sample_type,median2), color= Expression)) + 
  geom_boxplot() + 
  labs( x = "Expression (log2 normalized TPM)", y = "GTEx",title = "Isoform 2") + 
  theme(axis.text.y = element_text(size = 12), axis.title = element_text(size = 14), 
        plot.title = element_text( face = "bold", hjust = 0.5)) + 
  geom_vline(xintercept= c(iso1upper_bound, iso1upper_bound) , linetype="dashed", color = "red")


##############################################################################################################
################################################## Figure 2 ##################################################
######### Fig.2A #######
# Data visualization -- compare the isoform 1 and 2 expression in general between GTEX and TCGA data 
# Isoform 1 
merge1 %>% 
  ggplot(aes(Project,ENST00000620340.4, color = Project)) + geom_boxplot() + 
  stat_compare_means(method ="wilcox.test", label.y = 4, label.x = 1.25, size =4.5) + 
  theme_light() + my_theme +  theme(legend.position = "top", legend.background = element_rect(fill = "white")) + 
  scale_color_manual(labels = c("GTEX","TCGA"), values =c("#00BFC4","#F8766D"))

# Isoform 2
merge1 %>% 
  ggplot(aes(Project,ENST00000262752.4, color = Project)) + geom_boxplot() + 
  stat_compare_means(method ="wilcox.test", label.y = 3.5, label.x = 1.25, size =4.5) + 
  theme_light() + my_theme +  theme(legend.position = "top", legend.background = element_rect(fill = "white")) + 
  scale_color_manual(labels = c("GTEX","TCGA"), values =c("#00BFC4","#F8766D"))

merge1 %>% 
  ggplot(aes(Project,ENST00000262752.4, color = Project)) + geom_boxplot() +  
  theme_light() + my_theme + theme(legend.position = "top", legend.background = element_rect(fill = "white")) + 
  scale_color_manual(labels = c("GTEX","TCGA"), values =c("#00BFC4","#F8766D")) 

######### Fig.2B #######
## Pair-wise comparison of RSK4 expression between normal and tumour samples 
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
  geom_boxplot() + #ENST00000262752.4
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1)) + theme(plot.title = element_text(hjust = 0.5)) + facet_wrap(~ Group)


# Add p-value 
p + stat_compare_means(method ="wilcox.test", label.y = 1.5, label.x = 1.5, label = "p.signif") + 
  theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12),
        legend.text = element_text(size = 10), legend.title = element_text(size =12),
        strip.text = element_text(size = 14, face = "bold")) + theme_linedraw() + 
  scale_color_manual(labels = c("GTEX","TCGA"), values =c("#00BFC4","#F8766D")) + 
  theme(panel.grid.major = element_line(colour = "grey80"),panel.grid.minor = element_line(colour = "grey80"))


##############################################################################################################
################################################## Figure 3 ##################################################
######### Fig.3B #######
# Survival analysis  (Kaplan Meier curve)
## LGG - iso1
LGG.KM <- merge3 %>% filter(cancer.type.abbreviation == "LGG", !is.na(OS.time.x)) %>% 
  mutate(iso1.high = as.numeric(ENST00000620340.4 > median(ENST00000620340.4)),
         iso2.high = as.numeric(ENST00000262752.4 > median(ENST00000262752.4)))

LGG.os <- Surv(LGG.KM$OS.time.x, LGG.KM$OS.x)

ggsurvplot(survfit(LGG.os ~ LGG.KM$iso1.high), data = LGG.KM, conf.int = TRUE, 
           risk.table = TRUE, 
           risk.table.col = "strata",
           legend.title = "RSK4",
           legend.labs = c(paste0("Isoform 1-low"), paste0("Isoform 1-high")),
           palette = c("#00BFC4","#F8766D"),
           xlab = "OS time (days)",
           ylab = "Proportional survival rate",
           ggtheme = theme_bw() + my_theme)


## STAD - iso1 
STAD.KM <- merge3 %>% filter(cancer.type.abbreviation == "STAD", !is.na(OS.time.x), !is.na(OS.x)) %>% 
  mutate(iso1.high = as.numeric(ENST00000262752.4 > median(ENST00000262752.4)),
         iso2.high = as.numeric(ENST00000262752.4 > median(ENST00000262752.4)))

STAD.os <- Surv(STAD.KM$OS.time.x, STAD.KM$OS.x)

ggsurvplot(survfit(STAD.os ~ STAD.KM$iso1.high), data = STAD.KM, conf.int = TRUE, 
           risk.table = TRUE, 
           risk.table.col = "strata",
           legend.title = "RSK4",
           legend.labs = c(paste0("Isoform 1-low"), paste0("Isoform 1-high")),
           palette = c("#00BFC4","#F8766D"),
           xlab = "OS time (days)",
           ylab = "Proportional survival rate",
           ggtheme = theme_bw() + my_theme)

## CESC - iso2
CESC.KM <- merge3 %>% filter(cancer.type.abbreviation == "CESC", !is.na(OS.time.x), !is.na(OS.x)) %>% 
  mutate(iso1.high = as.numeric(ENST00000620340.4 > median(ENST00000620340.4)),
         iso2.high = as.numeric(ENST00000262752.4 > median(ENST00000262752.4)))

CESC.os <- Surv(CESC.KM$OS.time.x, CESC.KM$OS.x)

ggsurvplot(survfit(CESC.os ~ CESC.KM$iso2.high), data = CESC.KM, conf.int = TRUE, 
           risk.table = TRUE, 
           risk.table.col = "strata",
           legend.title = "RSK4",
           legend.labs = c(paste0("Isoform 1-low"), paste0("Isoform 1-high")),
           palette = c("#00BFC4","#F8766D"),
           xlab = "OS time (days)",
           ylab = "Proportional survival rate",
           ggtheme = theme_bw() + my_theme)

## KIRC - iso2
KIRC.KM <- merge3 %>% filter(cancer.type.abbreviation == "KIRC", !is.na(OS.time.x), !is.na(OS.x)) %>% 
  mutate(iso1.high = as.numeric(ENST00000620340.4 > median(ENST00000620340.4)),
         iso2.high = as.numeric(ENST00000262752.4 > median(ENST00000262752.4)))

KIRC.os <- Surv(KIRC.KM$OS.time.x, KIRC.KM$OS.x)

ggsurvplot(survfit(KIRC.os ~ KIRC.KM$iso2.high), data = KIRC.KM, conf.int = TRUE, 
           risk.table = TRUE, 
           risk.table.col = "strata",
           legend.title = "RSK4",
           legend.labs = c(paste0("Isoform 1-low"), paste0("Isoform 1-high")),
           palette = c("#00BFC4","#F8766D"),
           xlab = "OS time (days)",
           ylab = "Proportional survival rate",
           ggtheme = theme_bw() + my_theme)

##############################################################################################################
################################################## Figure 4 ##################################################
######### Fig.4A #######
# Correlation of isoform expression with clinical features
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

## LGG - iso 1 with cancer recurrence 
LGG.rsk <- merge_combine %>% 
  filter(Sample_type == "Brain Lower Grade Glioma") %>% 
  select(ENST00000620340.4,ENST00000262752.4, ratio, OS, OS.time)

LGG.merge <- merge(LGG.rsk,LGG.pheno,by = 0)

# Covariate 5: Cancer recurrence 
recur <- LGG.merge$new_tumor_event_after_initial_treatment
table(recur)
yes <- grep("YES",recur)
recur <- rep(0,nrow(LGG.merge))
recur[yes] <- 1

ggplot(data = LGG.merge, aes(as.factor(recur),ENST00000620340.4, color=as.factor(recur))) + geom_boxplot() + 
  stat_compare_means(method = "wilcox.test", paired = FALSE, label.x = 1.25) + # high stage show lower isoform 1 expression
  theme_bw() + 
  labs(x = "Cancer Recurrence", y = "Isoform 1 expression") + 
  scale_color_manual(values = c("blue","red")) +
  theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14),
        plot.title = element_text(size = 16), legend.position = "NONE") + my_theme

## KIRC - iso 2 with pathologic stage and metastasis
KIRC.rsk <- merge_combine %>% 
  filter(Sample_type == "Kidney Clear Cell Carcinoma") %>% 
  select(ENST00000620340.4,ENST00000262752.4, ratio, OS, OS.time)

KIRC.merge <- merge(KIRC.rsk,KIRC.pheno,by = 0)

# Covariate 3: Pathologic stage
stage <- KIRC.merge$pathologic_stage
table(stage)
x3<- grep("III",stage)
x4<- grep("IV",stage)
stage.high<-rep(0,nrow(KIRC.merge))
stage.high[c(x3,x4)]<-1

ggplot(KIRC.merge, aes(as.factor(stage.high),ENST00000262752.4, color = as.factor(stage.high) )) + geom_boxplot()+ 
  stat_compare_means(method = "wilcox.test", paired = FALSE, label.x = 1.25) + # high stage show lower isoform 2 expression
  theme_bw() + 
  labs(x = "Stage.high", y = "Isoform 2 expression") + 
  scale_color_manual(values = c("blue","red")) +
  theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14),
        plot.title = element_text(size = 16), legend.position = "NONE") + my_theme

# Covariate 5: Metastasis 
meta <- KIRC.merge$pathologic_M
table(meta)
M1 <- grep("M1",meta)
metastasis <- rep(0,nrow(KIRC.merge))
metastasis[M1] <- 1 

ggplot(KIRC.merge, aes(as.factor(metastasis),ENST00000262752.4, color = as.factor(metastasis))) + geom_boxplot()+ 
  stat_compare_means(method = "wilcox.test", paired = FALSE, label.x = 1.25) + # high stage show lower isoform 2 expression
  theme_bw() + 
  labs(x = "Metastasis", y = "Isoform 2 expression") + 
  scale_color_manual(values = c("blue","red")) +
  theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14),
        plot.title = element_text(size = 16), legend.position = "NONE") + my_theme

######### Fig.4B #######
## Barplot showing the relative proportion of immune cells in cancer 
# LGG 
LGG.merge <- merge(LGG.epic1, LGG.rsk, by = 0)
LGG.epic2 <- LGG.epic1[1:20,]
LGG.epic2 %>%
  rownames_to_column( var = "Samples") %>% 
  gather(Cell_type, fraction, B.cell:NK.cell) %>% 
  ggplot(aes(x = Samples, y = fraction, fill = Cell_type)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_y_continuous(labels = scales::percent_format()) + 
  labs(y = "Relative fraction", x = "Samples") + 
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1), 
        plot.title = element_text(hjust = 0.5))  # remove legends 

# STAD
STAD.merge <- merge(STAD.epic1, STAD.rsk, by = 0)
STAD.epic2 <- STAD.epic1[1:20,]
STAD.epic2 %>%
  rownames_to_column( var = "Samples") %>% 
  gather(Cell_type, fraction, B.cell:NK.cell) %>% 
  ggplot(aes(x = Samples, y = fraction, fill = Cell_type)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_y_continuous(labels = scales::percent_format()) + 
  labs(y = "Relative fraction", x = "Samples") + 
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1), 
        plot.title = element_text(hjust = 0.5))  # remove legends 

######### Fig.4C #######
# Survival plot demonstrate the association between immune cell fractions and patients' survival 
# LGG with CD4+ T cells 
LGG.os <- Surv(LGG.merge$OS.time, LGG.merge$OS)
summary(coxph(LGG.os ~ LGG.merge$T.cell.CD4.))
CD4.high <- as.numeric(LGG.merge$T.cell.CD4.>median(LGG.merge$T.cell.CD4.))

fit <- survfit(LGG.os ~ CD4.high)
legend.names <- c(paste0("CD4 T cell","-low"), paste0("CD4 T cell","-high")) # IF / ratio 
ggsurvplot(fit, data= LGG.os,
           pval = TRUE, pval.method = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           risk.table.col = "strata", # Change risk table color by group 
           #surv.median.line = "hv", # Specify median survival
           legend.title = "",
           legend.labs= legend.names,
           palette = c("#00BFC4","#F8766D"),
           xlab= "OS time (days)",
           ylab= "Proportion survival rate",
           ggtheme = theme_bw() + my_theme)

# STAD with CAFs 
STAD.os <- Surv(STAD.merge$OS.time, STAD.merge$OS)
summary(coxph(STAD.os ~ STAD.merge$CAFs))
CAFs.high <- as.numeric(STAD.merge$CAFs>median(STAD.merge$CAFs))

fit <- survfit(STAD.os ~ CAFs.high)
legend.names <- c(paste0("CAFs","-low"), paste0("CAFs","-high")) # IF / ratio 
ggsurvplot(fit, data= STAD.os,
           pval = TRUE, pval.method = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           risk.table.col = "strata", # Change risk table color by group 
           #surv.median.line = "hv", # Specify median survival
           legend.title="",
           legend.labs= legend.names,
           palette = c("#00BFC4","#F8766D"),
           xlab= "OS time (days)",
           ylab= "Proportion survival rate",
           ggtheme = theme_bw() + my_theme)

######### Fig.4D #######
# Correlation between isoform expression and relative proportion of immune cells 

# LGG with CD4+ T cells (isoform 1) 
ggplot(LGG.merge, aes(ENST00000620340.4, T.cell.CD4.)) + 
  geom_point() + stat_smooth(method = "lm") + stat_cor() + theme_bw() + my_theme + 
  labs(x = "Isoform 1 expression", y = "Proportion of CD4+ T cells")

# Non-linear regression 
ggplot(LGG.merge, aes(ENST00000620340.4, T.cell.CD4.)) + 
  geom_point() + stat_smooth(method = "nls", formula = 'y~a*exp(b*x)',method.args = list(start=c(a=0.1646, b=9.5e-8)), se=FALSE) + 
  stat_cor() + theme_bw() + my_theme + 
  labs(x = "Isoform 1 expression", y = "Proportion of CD4+ T cells")

# Take the top and bottom 25% of LGG patients by proportion of CD4+ T cells
LGG.merge %>% filter(T.cell.CD4. > quantile(T.cell.CD4., 0.75)) %>% count()

LGG.merge %>% filter(ENST00000620340.4 < quantile(ENST00000620340.4, 0.25)) %>% count()

LGG.merge$group <- NA
High <- which(LGG.merge$T.cell.CD4. > quantile(LGG.merge$T.cell.CD4., 0.75))
Low<- which(LGG.merge$ENST00000620340.4 < quantile(LGG.merge$ENST00000620340.4, 0.25))
overlap <- intersect(High, Low)
High <- High[-which(High %in% overlap)]
Low <- Low[-which(Low %in% overlap)]
LGG.merge$group[High] <- "high"
LGG.merge$group[Low] <- "low"

p <- ggplot(LGG.merge, aes(ENST00000620340.4, T.cell.CD4.)) + 
  geom_point() + stat_smooth(method = "nls", formula = 'y~a*exp(b*x)',method.args = list(start=c(a=0.1646, b=9.5e-8)), se=FALSE) + 
  stat_cor() + theme_bw() + my_theme + 
  labs(x = "Isoform 1 expression", y = "Proportion of CD4+ T cells")
p + geom_point(aes(ENST00000620340.4, T.cell.CD4., color = group)) 


fit <- survfit(STAD.os ~ CAFs.high)
legend.names <- c(paste0("CAFs","-low"), paste0("CAFs","-high")) # IF / ratio 
ggsurvplot(survfit(LGG.os ~ LGG.merge$group), data= LGG.merge,
           pval = TRUE, pval.method = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           risk.table.col = "strata", # Change risk table color by group 
           #surv.median.line = "hv", # Specify median survival
           legend.title="LGG",
           #legend.labs= legend.names,
           palette = c("#00BFC4","#F8766D"),
           xlab= "OS time (days)",
           ylab= "Proportion survival rate",
           ggtheme = theme_bw() + my_theme)

# STAD with CAFs 
ggplot(STAD.merge, aes(ENST00000620340.4, CAFs)) + 
  geom_point() + geom_smooth(method = "lm") + stat_cor() + theme_bw() + my_theme + 
  labs(x = "Isoform 1 expression", y = "Proportion of CAFs")
