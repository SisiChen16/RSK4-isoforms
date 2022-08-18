# Aim: To correlate immune cell infiltration with patients' overall survival and RSK4 isoform expression 

# Load packages 
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(dutchmasters)
library(broom)
library(RColorBrewer)
library(survival)
library(plotly)
library(viridis)
library(ggpubr)
library(survival)
library(survminer)

# Load data 
load("merge_combine.RData")

cancers <- c("Brain Lower Grade Glioma","Stomach Adenocarcinoma","Kidney Clear Cell Carcinoma","Cervical & Endocervical Cancer")

cancer.list <- c("LGG","STAD","KIRC","CESC")

ciber.list <- list(LGG.ciber1,STAD.ciber1,KIRC.ciber1,CESC.ciber1)
xcell.list <- list(LGG.xcell1,STAD.xcell1,KIRC.xcell1,CESC.xcell1)
epic.list <- list(LGG.epic1,STAD.epic1,KIRC.epic1,CESC.epic1)
quantiseq.list <- list(LGG.quantiseq1,STAD.quantiseq1,KIRC.quantiseq1,CESC.quantiseq1)
mcp.list <- list(LGG.mcp1,STAD.mcp1,KIRC.mcp1,CESC.mcp1)

rsk.list <- list(LGG.rsk, STAD.rsk, KIRC.rsk, CESC.rsk)


## 1. Correlate immune cell infiltration with overall survival 

ciber.surv.list <- list()
xcell.surv.list <- list()
epic.surv.list <- list()
quantiseq.surv.list <- list()
mcp.surv.list <- list()

for (i in 1:length(cancers)){
  ## 1. Prepare the data 
  # Load the RSK4 expression and clinical data for each cancer type 
  rsk.data <- merge_combine %>% 
    filter(Sample_type == cancers[i]) %>% 
    select(ENST00000620340.4,ENST00000262752.4,IF1, IF2, ratio, OS,OS.time)
  
  assign(paste0(cancer.list[i],".rsk"), data.frame(rsk.data))
  
  
  
  ## 2. Survival analysis -- by each deconvolution algorithm
  ### Cibersort
  ciber <- as.data.frame(ciber.list[[i]])
  merge <- merge(ciber, rsk.data, by = 0)
  merge <- merge %>% 
    column_to_rownames(var = "Row.names")
  
  ciber.result <- array(NA, c(6,4))
  colnames(ciber.result) <- c("HR","LCI","UCI","PVAL")
  rownames(ciber.result) <- colnames(merge)[1:6]
  ciber.result <- as.data.frame(ciber.result)
  
  # Create the survival objects
  os <- Surv(rsk.data$OS.time, rsk.data$OS)
  assign(paste0(cancer.list[i],".os"),os)
  
  for (j in 1:6){
    coxphmodel <- coxph(os ~ merge[,j])
    
    ciber.result$HR[j] <- summary(coxphmodel)$coef[1,2]
  
    ciber.result$LCI[j] <- summary(coxphmodel)$conf.int[1,3]
  
    ciber.result$UCI[j] <- summary(coxphmodel)$conf.int[1,4]
  
    ciber.result$PVAL[j] <- summary(coxphmodel)$coef[1,5]
  }
  
  ciber.result$Cancer <- cancer.list[i]
  ciber.result$method <- "CIBERSORT"
  ciber.surv.list[[i]] <- ciber.result
  # assign(paste0(cancer.list[i], ".ciber.surv"), data.frame(ciber.result))
  
  ### xCell
  xcell <- as.data.frame(xcell.list[[i]])
  xcell.merge <- merge(xcell, rsk.data, by = 0)
  xcell.merge <- xcell.merge %>% 
    column_to_rownames(var = "Row.names")
  xcell.result <- array(NA, c(6,4))
  colnames(xcell.result) <- c("HR","LCI","UCI","PVAL")
  rownames(xcell.result) <- colnames(xcell.merge)[1:6]
  xcell.result <- as.data.frame(xcell.result)
  
  for (j in 1:6){
    coxphmodel <- coxph(os ~ xcell.merge[,j])
    
    xcell.result$HR[j] <- summary(coxphmodel)$coef[1,2]
  
    xcell.result$LCI[j] <- summary(coxphmodel)$conf.int[1,3]
  
    xcell.result$UCI[j] <- summary(coxphmodel)$conf.int[1,4]
  
    xcell.result$PVAL[j] <- summary(coxphmodel)$coef[1,5]
  }
  
  xcell.result$Cancer <- cancer.list[i]
  xcell.result$method <- "xCell"
  xcell.surv.list[[i]] <- xcell.result
  
  ### EPIC
  epic <- as.data.frame(epic.list[[i]])
  epic.merge <- merge(epic, rsk.data, by = 0)
  epic.merge <- epic.merge %>% 
    column_to_rownames(var = "Row.names")
  epic.result <- array(NA, c(6,4))
  colnames(epic.result) <- c("HR","LCI","UCI","PVAL")
  rownames(epic.result) <- colnames(epic.merge)[1:6]
  epic.result <- as.data.frame(epic.result)
  
  for (j in 1:6){
    coxphmodel <- coxph(os ~ epic.merge[,j])
    
    epic.result$HR[j] <- summary(coxphmodel)$coef[1,2]
  
    epic.result$LCI[j] <- summary(coxphmodel)$conf.int[1,3]
  
    epic.result$UCI[j] <- summary(coxphmodel)$conf.int[1,4]
  
    epic.result$PVAL[j] <- summary(coxphmodel)$coef[1,5]
  }
  
  epic.result$Cancer <- cancer.list[i]
  epic.result$method <- "EPIC"
  epic.surv.list[[i]] <- epic.result
  
  ### QuanTIseq
  quantiseq <- as.data.frame(quantiseq.list[[i]])
  quantiseq.merge <- merge(quantiseq, rsk.data, by = 0)
  quantiseq.merge <- quantiseq.merge %>% 
    column_to_rownames(var = "Row.names")
  quantiseq.result <- array(NA, c(6,4))
  colnames(quantiseq.result) <- c("HR","LCI","UCI","PVAL")
  rownames(quantiseq.result) <- colnames(quantiseq.merge)[1:6]
  quantiseq.result <- as.data.frame(quantiseq.result)
  
  for (j in 1:6){
    coxphmodel <- coxph(os ~ quantiseq.merge[,j])
    
    quantiseq.result$HR[j] <- summary(coxphmodel)$coef[1,2]
  
    quantiseq.result$LCI[j] <- summary(coxphmodel)$conf.int[1,3]
  
    quantiseq.result$UCI[j] <- summary(coxphmodel)$conf.int[1,4]
  
    quantiseq.result$PVAL[j] <- summary(coxphmodel)$coef[1,5]
  }
  
  quantiseq.result$Cancer <- cancer.list[i]
  quantiseq.result$method <- "quanTIseq"
  quantiseq.surv.list[[i]] <- quantiseq.result
  
  
  ### mcpCounter
  mcp <- as.data.frame(mcp.list[[i]])
  mcp.merge <- merge(mcp, rsk.data, by = 0)
  mcp.merge <- mcp.merge %>% 
    column_to_rownames(var = "Row.names")
  mcp.result <- array(NA, c(6,4))
  colnames(mcp.result) <- c("HR","LCI","UCI","PVAL")
  rownames(mcp.result) <- colnames(mcp.merge)[1:6]
  mcp.result <- as.data.frame(mcp.result)
  
  for (j in 1:6){
    coxphmodel <- coxph(os ~ mcp.merge[,j])
    
    mcp.result$HR[j] <- summary(coxphmodel)$coef[1,2]
  
    mcp.result$LCI[j] <- summary(coxphmodel)$conf.int[1,3]
  
    mcp.result$UCI[j] <- summary(coxphmodel)$conf.int[1,4]
  
    mcp.result$PVAL[j] <- summary(coxphmodel)$coef[1,5]
  }

  mcp.result$Cancer <- cancer.list[i]
  mcp.result$method <- "mcpCounter"
  mcp.surv.list[[i]] <- mcp.result
}


ciber.surv.list
xcell.surv.list
epic.surv.list
quantiseq.surv.list
mcp.surv.list

save(ciber.surv.list, file = "ciber.surv.list.RData")
save(xcell.surv.list, file = "xcell.surv.list.RData")
save(epic.surv.list, file = "epic.surv.list.RData")
save(quantiseq.surv.list, file = "quantiseq.surv.list.RData")
save(mcp.surv.list, file = "mcp.surv.list.RData")


## 2. Correlate immune cell infiltration with RSK4 isoform expression 

#### Cibersort 
ciber.rsk.list <- list()
for (i in 1:length(cancers)){
    # Import and merge dataframe
  rsk.data <- as.data.frame(rsk.list[[i]])
  ciber <- as.data.frame(ciber.list[[i]])
  merge <- merge(ciber,rsk.data, by = 0)
  merge <- merge %>% 
    column_to_rownames(var = "Row.names")
  
    # Correlation test (isoform 1, 2 and ratio)
  iso1.result <- data.frame()
  iso2.result <- data.frame()
  ratio.result <- data.frame()
  for (j in 1:6){
    # Isoform 1 correlation 
    iso1.cor <- tidy(cor.test(merge$ENST00000620340.4, merge[,j],method = "spearman",
                          continuity = TRUE, exact = FALSE))
    iso1.cor$cell_type <- colnames(merge)[j]
    
    # Isoform 2 correlation 
    iso2.cor <- tidy(cor.test(merge$ENST00000262752.4, merge[,j],method = "spearman",
                          continuity = TRUE, exact = FALSE))
    iso2.cor$cell_type <- colnames(merge)[j]
    
    # ratio correlation 
    ratio.cor <- tidy(cor.test(merge$ratio, merge[,j],method = "spearman",
                          continuity = TRUE, exact = FALSE))
    ratio.cor$cell_type <- colnames(merge)[j]

    iso1.result <- rbind(iso1.result, iso1.cor)
    iso2.result <- rbind(iso2.result, iso2.cor)
    ratio.result <- rbind(ratio.result, ratio.cor)
  }
  iso1.result$parameter <- "Isoform 1"
  iso2.result$parameter <- "Isoform 2"
  ratio.result$parameter <- "ratio"
  
  merge.result <- rbind(iso1.result, iso2.result, ratio.result)
  merge.result$cancer_type <- cancer.list[i]
  merge.result$cell_decon <- "CIBERSORT"
  ciber.rsk.list[[i]] <- merge.result
  
}

for (i in 1:length(ciber.rsk.list)){
  results <- as.data.frame(ciber.rsk.list[[i]])
  results <- results %>% 
    filter(p.value < 0.05 )
  name <- unique(results$cancer_type)
  
  assign(paste0(name,".ciber.rsk"), data.frame(results))
}

LGG.ciber.rsk
STAD.ciber.rsk 
KIRC.ciber.rsk
CESC.ciber.rsk

#### xCell 
xcell.rsk.list <- list()

for (i in 1:length(cancers)){
    # Import and merge dataframe
  rsk.data <- as.data.frame(rsk.list[[i]])
  xcell <- as.data.frame(xcell.list[[i]])
  merge.xcell <- merge(xcell,rsk.data, by = 0)
  merge.xcell <- merge.xcell %>% 
    column_to_rownames(var = "Row.names")
  
    # Correlation test (isoform 1, 2 and ratio)
  iso1.result <- data.frame()
  iso2.result <- data.frame()
  ratio.result <- data.frame()
  for (j in 1:6){
    # Isoform 1 correlation 
    iso1.cor <- tidy(cor.test(merge.xcell$ENST00000620340.4, merge.xcell[,j],method = "spearman",
                          continuity = TRUE, exact = FALSE))
    iso1.cor$cell_type <- colnames(merge.xcell)[j]
    
    # Isoform 2 correlation 
    iso2.cor <- tidy(cor.test(merge.xcell$ENST00000262752.4, merge.xcell[,j],method = "spearman",
                          continuity = TRUE, exact = FALSE))
    iso2.cor$cell_type <- colnames(merge.xcell)[j]
    
    # ratio correlation 
    ratio.cor <- tidy(cor.test(merge.xcell$ratio, merge.xcell[,j],method = "spearman",
                          continuity = TRUE, exact = FALSE))
    ratio.cor$cell_type <- colnames(merge.xcell)[j]

    iso1.result <- rbind(iso1.result, iso1.cor)
    iso2.result <- rbind(iso2.result, iso2.cor)
    ratio.result <- rbind(ratio.result, ratio.cor)
  }
  iso1.result$parameter <- "Isoform 1"
  iso2.result$parameter <- "Isoform 2"
  ratio.result$parameter <- "ratio"
  
  merge.result <- rbind(iso1.result, iso2.result, ratio.result)
  merge.result$cancer_type <- cancer.list[i]
  merge.result$cell_decon <- "xCell"
  xcell.rsk.list[[i]] <- merge.result
  
}


for (i in 1:length(xcell.rsk.list)){
  results <- as.data.frame(xcell.rsk.list[[i]])
  results <- results %>% 
    filter(p.value < 0.05 )
  name <- unique(results$cancer_type)
  
  assign(paste0(name,".xcell.rsk"), data.frame(results))
}


LGG.xcell.rsk
STAD.xcell.rsk 
KIRC.xcell.rsk
CESC.xcell.rsk 


#### EPIC 
epic.rsk.list <- list()
for (i in 1:length(cancers)){
    # Import and merge dataframe
  rsk.data <- as.data.frame(rsk.list[[i]])
  epic <- as.data.frame(epic.list[[i]])
  merge.epic <- merge(epic,rsk.data, by = 0)
  merge.epic <- merge.epic %>% 
    column_to_rownames(var = "Row.names")
  
    # Correlation test (isoform 1, 2 and ratio)
  iso1.result <- data.frame()
  iso2.result <- data.frame()
  ratio.result <- data.frame()
  for (j in 1:6){
    # Isoform 1 correlation 
    iso1.cor <- tidy(cor.test(merge.epic$ENST00000620340.4, merge.epic[,j],method = "spearman",
                          continuity = TRUE, exact = FALSE))
    iso1.cor$cell_type <- colnames(merge.epic)[j]
    
    # Isoform 2 correlation 
    iso2.cor <- tidy(cor.test(merge.epic$ENST00000262752.4, merge.epic[,j],method = "spearman",
                          continuity = TRUE, exact = FALSE))
    iso2.cor$cell_type <- colnames(merge.epic)[j]
    
    # ratio correlation 
    ratio.cor <- tidy(cor.test(merge.epic$ratio, merge.epic[,j],method = "spearman",
                          continuity = TRUE, exact = FALSE))
    ratio.cor$cell_type <- colnames(merge.epic)[j]

    iso1.result <- rbind(iso1.result, iso1.cor)
    iso2.result <- rbind(iso2.result, iso2.cor)
    ratio.result <- rbind(ratio.result, ratio.cor)
  }
  iso1.result$parameter <- "Isoform 1"
  iso2.result$parameter <- "Isoform 2"
  ratio.result$parameter <- "ratio"
  
  merge.result <- rbind(iso1.result, iso2.result, ratio.result)
  merge.result$cancer_type <- cancer.list[i]
  merge.result$cell_decon <- "EPIC"
  epic.rsk.list[[i]] <- merge.result
  
}


for (i in 1:length(epic.rsk.list)){
  results <- as.data.frame(epic.rsk.list[[i]])
  results <- results %>% 
    filter(p.value < 0.05 )
  name <- unique(results$cancer_type)
  
  assign(paste0(name,".epic.rsk"), data.frame(results))
}

LGG.epic.rsk
STAD.epic.rsk 
KIRC.epic.rsk
CESC.epic.rsk 


#### quanTIseq 
quantiseq.rsk.list <- list()

for (i in 1:length(cancers)){
    # Import and merge dataframe
  rsk.data <- as.data.frame(rsk.list[[i]])
  quantiseq <- as.data.frame(quantiseq.list[[i]])
  merge.quantiseq <- merge(quantiseq,rsk.data, by = 0)
  merge.quantiseq <- merge.quantiseq %>% 
    column_to_rownames(var = "Row.names")
  
    # Correlation test (isoform 1, 2 and ratio)
  iso1.result <- data.frame()
  iso2.result <- data.frame()
  ratio.result <- data.frame()
  for (j in 1:6){
    # Isoform 1 correlation 
    iso1.cor <- tidy(cor.test(merge.quantiseq$ENST00000620340.4, merge.quantiseq[,j],method = "spearman",
                          continuity = TRUE, exact = FALSE))
    iso1.cor$cell_type <- colnames(merge.quantiseq)[j]
    
    # Isoform 2 correlation 
    iso2.cor <- tidy(cor.test(merge.quantiseq$ENST00000262752.4, merge.quantiseq[,j],method = "spearman",
                          continuity = TRUE, exact = FALSE))
    iso2.cor$cell_type <- colnames(merge.quantiseq)[j]
    
    # ratio correlation 
    ratio.cor <- tidy(cor.test(merge.quantiseq$ratio, merge.quantiseq[,j],method = "spearman",
                          continuity = TRUE, exact = FALSE))
    ratio.cor$cell_type <- colnames(merge.quantiseq)[j]

    iso1.result <- rbind(iso1.result, iso1.cor)
    iso2.result <- rbind(iso2.result, iso2.cor)
    ratio.result <- rbind(ratio.result, ratio.cor)
  }
  iso1.result$parameter <- "Isoform 1"
  iso2.result$parameter <- "Isoform 2"
  ratio.result$parameter <- "ratio"
  
  merge.result <- rbind(iso1.result, iso2.result, ratio.result)
  merge.result$cancer_type <- cancer.list[i]
  merge.result$cell_decon <- "QuanTIseq"
  quantiseq.rsk.list[[i]] <- merge.result
  
}


for (i in 1:length(quantiseq.rsk.list)){
  results <- as.data.frame(quantiseq.rsk.list[[i]])
  results <- results %>% 
    filter(p.value < 0.05 )
  name <- unique(results$cancer_type)
  
  assign(paste0(name,".quantiseq.rsk"), data.frame(results))
}


LGG.quantiseq.rsk
STAD.quantiseq.rsk 
KIRC.quantiseq.rsk
CESC.quantiseq.rsk 




#### mcpCounter 

mcp.rsk.list <- list()

for (i in 1:length(cancers)){
    # Import and merge dataframe
  rsk.data <- as.data.frame(rsk.list[[i]])
  mcp <- as.data.frame(mcp.list[[i]])
  merge.mcp <- merge(mcp,rsk.data, by = 0)
  merge.mcp <- merge.mcp %>% 
    column_to_rownames(var = "Row.names")
  
    # Correlation test (isoform 1, 2 and ratio)
  iso1.result <- data.frame()
  iso2.result <- data.frame()
  ratio.result <- data.frame()
  for (j in 1:6){
    # Isoform 1 correlation 
    iso1.cor <- tidy(cor.test(merge.mcp$ENST00000620340.4, merge.mcp[,j],method = "spearman",
                          continuity = TRUE, exact = FALSE))
    iso1.cor$cell_type <- colnames(merge.mcp)[j]
    
    # Isoform 2 correlation 
    iso2.cor <- tidy(cor.test(merge.mcp$ENST00000262752.4, merge.mcp[,j],method = "spearman",
                          continuity = TRUE, exact = FALSE))
    iso2.cor$cell_type <- colnames(merge.mcp)[j]
    
    # ratio correlation 
    ratio.cor <- tidy(cor.test(merge.mcp$ratio, merge.mcp[,j],method = "spearman",
                          continuity = TRUE, exact = FALSE))
    ratio.cor$cell_type <- colnames(merge.mcp)[j]

    iso1.result <- rbind(iso1.result, iso1.cor)
    iso2.result <- rbind(iso2.result, iso2.cor)
    ratio.result <- rbind(ratio.result, ratio.cor)
  }
  iso1.result$parameter <- "Isoform 1"
  iso2.result$parameter <- "Isoform 2"
  ratio.result$parameter <- "ratio"
  
  merge.result <- rbind(iso1.result, iso2.result, ratio.result)
  merge.result$cancer_type <- cancer.list[i]
  merge.result$cell_decon <- "MCPcounter"
  mcp.rsk.list[[i]] <- merge.result
  
}


for (i in 1:length(mcp.rsk.list)){
  results <- as.data.frame(mcp.rsk.list[[i]])
  results <- results %>% 
    filter(p.value < 0.05 )
  name <- unique(results$cancer_type)
  
  assign(paste0(name,".mcp.rsk"), data.frame(results))
}

LGG.mcp.rsk
STAD.mcp.rsk 
KIRC.mcp.rsk
CESC.mcp.rsk 


######################## Data visualization ##########################
## Brain Lower Grade Glioma (LGG) 
LGG.merge <- merge(LGG.epic1, LGG.rsk, by= 0)

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

LGG.os <- Surv(LGG.merge$OS.time, LGG.merge$OS)

# Survival analysis
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



# CD4+ T cells (isoform 1,2 ratio)
ggplot(LGG.merge, aes(ENST00000620340.4, T.cell.CD4., color = OS.time)) + geom_point() + geom_smooth(method = "lm") + scale_color_viridis(option = "D") # ENST00000262752.4, ratio 

ggplot(LGG.merge, aes(ENST00000620340.4, T.cell.CD4.)) + 
  geom_point() + stat_smooth(method = "lm") + stat_cor() + theme_bw() + my_theme + 
  labs(x = "Isoform 1 expression", y = "Proportion of CD4+ T cells")

# Non-linear regression 
ggplot(LGG.merge, aes(ENST00000620340.4, T.cell.CD4.)) + 
  geom_point() + stat_smooth(method = "nls", formula = 'y~a*exp(b*x)',method.args = list(start=c(a=0.1646, b=9.5e-8)), se=FALSE) + 
  stat_cor() + theme_bw() + my_theme + 
  labs(x = "Isoform 1 expression", y = "Proportion of CD4+ T cells")

# Take the top 25% and bottom 25% of LGG patients for survival analysis
quantile(LGG.merge$T.cell.CD4., 0.75)
quantile(LGG.merge$ENST00000620340.4, 0.75)

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


summary(lm(LGG.merge$ENST00000620340.4 ~ LGG.merge$T.cell.CD4.))

ggplot(LGG.merge, aes(ENST00000262752.4, T.cell.CD4., color = OS.time)) + geom_point() + geom_smooth(method = "lm") + scale_color_viridis(option = "D")
summary(lm(LGG.merge$ENST00000262752.4 ~ LGG.merge$T.cell.CD4.))

ggplot(LGG.merge, aes(ratio, T.cell.CD4., color = OS.time)) + geom_point() + geom_smooth(method = "lm") + scale_color_viridis(option = "D")
summary(lm(LGG.merge$ratio ~ LGG.merge$T.cell.CD4.))

# CAFs
ggplot(LGG.merge, aes(ENST00000620340.4, CAFs, color = OS.time)) + geom_point() + geom_smooth(method = "lm") + scale_color_viridis(option = "D")
summary(lm(LGG.merge$ENST00000620340.4 ~ LGG.merge$CAFs))


### Adrenocortical Carcinoma (ACC)
ACC.merge <- merge(ACC.epic1, ACC.rsk, by = 0)
ACC.epic2 <- ACC.epic1[1:20,]
ACC.epic2 %>%
  rownames_to_column( var = "Samples") %>% 
  gather(Cell_type, fraction, B.cell:NK.cell) %>% 
  ggplot(aes(x = Samples, y = fraction, fill = Cell_type)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_y_continuous(labels = scales::percent_format()) + 
  labs(y = "Relative fraction", x = "Samples") + 
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1), 
        plot.title = element_text(hjust = 0.5))  # remove legends 

# CD4+ T cells (isoform 1,2 ratio)
ggplot(ACC.merge, aes(ratio, T.cell.CD4., color = OS.time)) + geom_point() + geom_smooth(method = "lm") + scale_color_viridis(option = "D")
summary(lm(ACC.merge$ratio ~ ACC.merge$T.cell.CD4.))


### Stomach adenocarcinoma (STAD)
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

# CD4+ T cells (isoform 1,2 ratio)
ggplot(STAD.merge, aes(ENST00000620340.4, T.cell.CD4., color = OS.time)) + geom_point() + geom_smooth(method = "lm") + scale_color_viridis(option = "D")
summary(lm(STAD.merge$ENST00000620340.4 ~ STAD.merge$T.cell.CD4.))

ggplot(STAD.merge, aes(ENST00000262752.4, T.cell.CD4., color = OS.time)) + geom_point() + geom_smooth(method = "lm") + scale_color_viridis(option = "D")
summary(lm(STAD.merge$ENST00000262752.4 ~ STAD.merge$T.cell.CD4.))

# CAFs 
ggplot(STAD.merge, aes(ENST00000620340.4, CAFs, color = OS.time)) + geom_point() + geom_smooth(method = "lm") + scale_color_viridis(option = "D")
summary(lm(STAD.merge$ENST00000620340.4 ~ STAD.merge$CAFs))

ggplot(STAD.merge, aes(ENST00000620340.4, CAFs)) + 
  geom_point() + geom_smooth(method = "lm") + stat_cor() + theme_bw() + my_theme + 
  labs(x = "Isoform 1 expression", y = "Proportion of CAFs")

ggplot(STAD.merge, aes(ENST00000262752.4, CAFs, color = OS.time)) + geom_point() + geom_smooth(method = "lm") + scale_color_viridis(option = "D")
summary(lm(STAD.merge$ENST00000262752.4 ~ STAD.merge$CAFs))

# Survival analysis
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


### Kidney Clear Cell Carcinoma (KIRC)
KIRC.merge <- merge(KIRC.epic1, KIRC.rsk, by = 0)
KIRC.epic2 <- KIRC.epic1[1:20,]
KIRC.epic2 %>%
  rownames_to_column( var = "Samples") %>% 
  gather(Cell_type, fraction, B.cell:NK.cell) %>% 
  ggplot(aes(x = Samples, y = fraction, fill = Cell_type)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_y_continuous(labels = scales::percent_format()) + 
  labs(y = "Relative fraction", x = "Samples") + 
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1), 
        plot.title = element_text(hjust = 0.5))  # remove legends 

# CD8+ T cells (isoform 2)
ggplot(KIRC.merge, aes(ENST00000262752.4, T.cell.CD8., color = OS.time)) + geom_point() + geom_smooth(method = "lm") + scale_color_viridis(option = "D")
summary(lm(KIRC.merge$ENST00000262752.4 ~ KIRC.merge$T.cell.CD8.))


# CD4+ T cells (isoform 1,2)
ggplot(KIRC.merge, aes(ENST00000620340.4, T.cell.CD4., color = OS.time)) + geom_point() + geom_smooth(method = "lm") + scale_color_viridis(option = "D")
summary(lm(KIRC.merge$ENST00000620340.4 ~ KIRC.merge$T.cell.CD4.))

ggplot(KIRC.merge, aes(ENST00000262752.4, T.cell.CD4., color = OS.time)) + geom_point() + geom_smooth(method = "lm") + scale_color_viridis(option = "D")
summary(lm(KIRC.merge$ENST00000262752.4 ~ KIRC.merge$T.cell.CD4.))


### Cervical & Endocervical carcinoma CESC 
CESC.merge <- merge(CESC.epic1, CESC.rsk, by = 0)
CESC.epic2 <- CESC.epic1[1:20,]
CESC.epic2 %>%
  rownames_to_column( var = "Samples") %>% 
  gather(Cell_type, fraction, B.cell:NK.cell) %>% 
  ggplot(aes(x = Samples, y = fraction, fill = Cell_type)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_y_continuous(labels = scales::percent_format()) + 
  labs(y = "Relative fraction", x = "Samples") + 
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1), 
        plot.title = element_text(hjust = 0.5))  # remove legends 


  
# CAFs
ggplot(CESC.merge, aes(ENST00000620340.4, CAFs, color = OS.time)) + geom_point() + geom_smooth(method = "lm") + scale_color_viridis(option = "D")
summary(lm(CESC.merge$ENST00000620340.4 ~ CESC.merge$CAFs))

ggplot(CESC.merge, aes(ENST00000262752.4, CAFs, color = OS.time)) + geom_point() + geom_smooth(method = "lm") + scale_color_viridis(option = "D")
summary(lm(CESC.merge$ENST00000262752.4 ~ CESC.merge$CAFs))

ggplot(CESC.merge, aes(log10(ratio), CAFs, color = OS.time)) + geom_point() + geom_smooth(method = "lm") + scale_color_viridis(option = "D")
summary(lm(log10(CESC.merge$ratio) ~ CESC.merge$CAFs))
