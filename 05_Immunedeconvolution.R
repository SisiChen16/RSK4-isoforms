# Aim: 
# 1. To use cell deconvolution algorithms to determine the immune cell landscape in each cancer samples
# 2. Standardize the outputs from five different algorithms by extracting the common cell types and calculate the relative cell fractions 

# Load packages 
library(dplyr)
library(ggplot2)
library(tidyr)
library(immunedeconv)
library(tibble)
library(dutchmasters)
library(broom)
library(RColorBrewer)
library(survival)
library(plotly)

# Input data 
setwd("/Users/chensisi/Documents/RNAseq/4_Cell_deconvolution")

LGG <- read.table("LGGexpression.tpm.txt", header = TRUE, row.names = 1, sep = "\t")
STAD <- read.table("STADexpression.tpm.txt", header = TRUE, row.names = 1, sep = "\t")
KIRC <- read.table("KIRCexpression.tpm.txt", header = TRUE, row.names = 1, sep = "\t")
CESC <- read.table("CESCexpression.tpm.txt", header = TRUE, row.names = 1, sep = "\t")


	# Convert the dataframe into matrix 
LGG <- as.matrix(LGG)
STAD <- as.matrix(STAD)
KIRC <- as.matrix(KIRC)
CESC <- as.matrix(CESC)


#save.image(file = "expression.RData")
load("expression.RData")

## Cibersort deconvolution 

source('Cibersort2.R')

LGG.ciber <- CIBERSORT('PanelA_signature_matrix.txt',"LGGexpression.tpm.txt", perm = 100, QN = T)

STAD.ciber <- CIBERSORT('PanelA_signature_matrix.txt',"STADexpression.tpm.txt", perm = 100, QN = T)

KIRC.ciber <- CIBERSORT('PanelA_signature_matrix.txt',"KIRCexpression.tpm.txt", perm = 100, QN = T)

CESC.ciber <- CIBERSORT('PanelA_signature_matrix.txt',"CESCexpression.tpm.txt", perm = 100, QN = T)


## EPIC deconvolution 
LGG.epic <- immunedeconv::deconvolute(LGG, "epic", tumor= TRUE )
STAD.epic <- immunedeconv::deconvolute(STAD,"epic", tumor= TRUE)
KIRC.epic <- immunedeconv::deconvolute(KIRC, "epic", tumor= TRUE)
CESC.epic <- immunedeconv::deconvolute(CESC, "epic", tumor= TRUE)


## xCell deconvolution 
LGG.xcell <- immunedeconv::deconvolute_xcell(LGG, array = FALSE) 
STAD.xcell <- immunedeconv::deconvolute_xcell(STAD, array = FALSE) 
KIRC.xcell <- immunedeconv::deconvolute_xcell(KIRC, array = FALSE) 
CESC.xcell <- immunedeconv::deconvolute_xcell(CESC, array = FALSE) 

## quanTIseq deconvolution 
LGG.quantiseq <- immunedeconv::deconvolute(LGG, "quantiseq", arrays = FALSE)
STAD.quantiseq <- immunedeconv::deconvolute(STAD, "quantiseq", arrays = FALSE)
KIRC.quantiseq <- immunedeconv::deconvolute(KIRC, "quantiseq", arrays = FALSE)
CESC.quantiseq <- immunedeconv::deconvolute(CESC, "quantiseq", arrays = FALSE)


## MCPcounter deconvolution 
LGG.mcp <- immunedeconv::deconvolute(LGG, "mcp_counter")
STAD.mcp <- immunedeconv::deconvolute(STAD, "mcp_counter")
KIRC.mcp <- immunedeconv::deconvolute(KIRC, "mcp_counter")
CESC.mcp <- immunedeconv::deconvolute(CESC, "mcp_counter")

## Use "For" loop to complete all cell deconvolution analysis 
inputs <- list(LGG, STAD, KIRC,CESC)
filenames <- c("LGGexpression.tpm.txt","STADexpression.tpm.txt","KIRCexpression.tpm.txt","CESCexpression.tpm.txt")
cancer_list <- c("LGG","STAD","KIRC","CESC")

for (i in 1:length(cancer_list)){
  
  print(filenames[i])
  print(cancer_list[i])
  data <- as.matrix(inputs[[i]])
  
  # Cibersort deconvolution 
  Cibersort <- CIBERSORT('PanelA_signature_matrix.txt',filenames[i], perm = 100, QN = T)
  
  # EPIC deconvolution
  Epic <- immunedeconv::deconvolute(data, "epic")
  assign(paste0(cancer.list[i],".epic"), data.frame(Epic))
  
  # xCell deconvolution 
  xCell <- immunedeconv::deconvolute_xcell(data, array = FALSE) 
  assign(paste0(cancer.list[i],".xcell"), data.frame(xCell))
  
  # Quantiseq deconvolution
  quantiseq <- immunedeconv::deconvolute(data, "quantiseq", arrays = FALSE)
  assign(paste0(cancer.list[i],".quantiseq"), data.frame(quantiseq))
  
  # MCPcounter deconvolution
  mcp <- immunedeconv::deconvolute(data, "mcp_counter")
  assign(paste0(cancer.list[i],".mcp"), data.frame(mcp))
  
}


## Standardize the outputs from all five algorithms by extracting the cell types shared in common and computing the relative cell fractions 

# Cibersort results 
inputs <- list(LGG.ciber, STAD.ciber, KIRC.ciber,CESC.ciber)
cancer.list <- c("LGG","STAD","KIRC","CESC")

for (i in 1:length(cancer.list)){
  # Discard irrelevant columns 
  discard <- c("P-value","Correlation","RMSE") 
  ciber.result <- as.data.frame(inputs[[i]])
  ciber.result1 <- ciber.result[,-which(colnames(ciber.result) %in% discard)]
  
  # Select the cell fractions of target cell type 
  ciber.result2 <- ciber.result1 %>%  # merge those cell types of interest 
    mutate(`B cell` = `naive B-cells`,
         `T cell CD4+` = `CD4+ T-cells`, 
         `T cell CD8+` = `CD8+ T-cells`, 
         `CAFs` = `Fibroblasts`, 
         `Macrophage` = `Macrophages` + `Macrophages M1` + `Macrophages M2`, 
         `NK cell` = `NK cells`) %>% 
  select("B cell","T cell CD4+", "T cell CD8+", "CAFs", "Macrophage","NK cell")

  # Export the results as data.frame 
  assign(paste0(cancer.list[i],".ciber1"), data.frame(ciber.result2))
}


# xCell results 

inputs <- list(KICH.xcell, LGG.xcell, ACC.xcell, STAD.xcell, KIRC.xcell,CESC.xcell)
cancer.list <- c("KICH","LGG","ACC","STAD","KIRC","CESC")


for (i in 1:length(inputs)){
  # Discard irrelevant columns from xcell results (immunescore, stromascore, microenvironmentScore)
  discard <- c("ImmuneScore","StromaScore","MicroenvironmentScore")
  xcell.result <- as.data.frame(t(inputs[[i]]))
  xcell.result1 <- xcell.result[,-which(colnames(xcell.result) %in% discard)]
  
  # Calculate the total value of enrichment scores from all cell types & calculate the fraction  
  xcell.result1$total <- rowSums(xcell.result1)
  xcell.result2 <- xcell.result1 %>%  # merge those cell types of interest 
    mutate(`B cell` = `B-cells` + `Memory B-cells` + `naive B-cells` + `pro B-cells`,
         `T cell CD4+` = `CD4+ memory T-cells` + `CD4+ naive T-cells` + `CD4+ T-cells` + `CD4+ Tcm` + `CD4+ Tem` + `Th1 cells` + `Tgd cells` + `Th2 cells` + `Tregs`, 
         `T cell CD8+` = `CD8+ naive T-cells` + `CD8+ Tcm` + `CD8+ T-cells` + `CD8+ Tem`, 
         `CAFs` = `Fibroblasts`, 
         `Macrophage` = `Macrophages` + `Macrophages M1` + `Macrophages M2`, 
         `NK cell` = `NK cells`) 
   
  # Calculate the fraction of each cell type
  xcell.result3 <- xcell.result2 %>% 
    mutate(`B cell` = `B cell`/total, 
         `T cell CD4+` = `T cell CD4+`/total, 
         `T cell CD8+` = `T cell CD8+`/total,
         `CAFs` = `CAFs`/total,
         `Macrophage` = `Macrophage`/total,
         `NK cell` = `NK cell`/total)
  
  # Select the cell types of interest 
  cell_types <- c("B cell","T cell CD4+", "T cell CD8+", "CAFs", "Macrophage","NK cell")
  xcell.result4 <- xcell.result3[, which(colnames(xcell.result3) %in% cell_types)]

  assign(paste0(cancer.list[i],".xcell1"), data.frame(xcell.result4))
}


# EPIC results 

inputs <- list(KICH.epic, LGG.epic, ACC.epic, STAD.epic, KIRC.epic,CESC.epic)
cancer.list <- c("KICH","LGG","ACC","STAD","KIRC","CESC")

for (i in 1:length(cancer.list)){
  # Input data and convert column to rownames 
  epic.result <- as.data.frame(inputs[[i]])
  epic.result1 <- epic.result %>% 
    column_to_rownames(var = "cell_type")
  
  # Transpose the dataframe 
  epic.result2 <- t(epic.result1)
  rowSums(epic.result2) # Equals to 1
  
  # Select cell types of interest 
  epic.result3 <- as.data.frame(epic.result2) %>% 
    mutate(`CAFs` = `Cancer associated fibroblast`) %>% 
    select("B cell","T cell CD4+", "T cell CD8+", "CAFs", "Macrophage","NK cell")
  
  # Export the results as data.frame 
  assign(paste0(cancer.list[i],".epic1"), data.frame(epic.result3))
}


# quanTIseq results 
inputs <- list(KICH.quantiseq, LGG.quantiseq, ACC.quantiseq, STAD.quantiseq, KIRC.quantiseq,CESC.quantiseq)
cancer.list <- c("KICH","LGG","ACC","STAD","KIRC","CESC")

for (i in 1:length(cancer.list)){
  # Input data and convert column to rownames 
  quantiseq.result <- as.data.frame(inputs[[i]])
  quantiseq.result1 <- quantiseq.result %>% 
    column_to_rownames(var = "cell_type")
  
  # Transpose the dataframe 
  quantiseq.result2 <- t(quantiseq.result1)
  rowSums(quantiseq.result2) # Equals to 1 -- make sure the cell fractions add up to 1 
  
  # Select the cell types of interest 
  quantiseq.result3 <- as.data.frame(quantiseq.result2) %>% 
    mutate(`Macrophage` = `Macrophage M1` + `Macrophage M2`,
           `T cell CD4+` = `T cell CD4+ (non-regulatory)` + `T cell regulatory (Tregs)`) %>% 
    select("B cell","T cell CD4+", "T cell CD8+", "Macrophage","NK cell") # CAFs are lacking  
  
  # Export the results as data.frame 
  assign(paste0(cancer.list[i],".quantiseq1"), data.frame(quantiseq.result3))
  
}


# mcpCounter results 

inputs <- list(KICH.mcp, LGG.mcp, ACC.mcp, STAD.mcp, KIRC.mcp,CESC.mcp)
cancer.list <- c("KICH","LGG","ACC","STAD","KIRC","CESC")

for (i in 1:length(cancer.list)){
  # Input data and convert column to rownames 
  mcp.result <- as.data.frame(inputs[[i]])
  mcp.result1 <- mcp.result %>% 
    column_to_rownames(var = "cell_type")
  mcp.result1 <- t(mcp.result1)
  
  # Discard irrelevant columns 
  discard <- c("cytotoxicity score") 
  mcp.result2 <- as.data.frame(mcp.result1[,-which(colnames(mcp.result1) %in% discard)])
  
  # Select the cell types of interest & Calculate the cell fractions 
  colnames(mcp.result2)
  mcp.result2$total <- rowSums(mcp.result2)
  mcp.result3 <- as.data.frame(mcp.result2) %>% 
    mutate(`B cell` = `B cell`/total,
           `T cell CD4+` = `T cell`/total, 
           `T cell CD8+` = `T cell CD8+`/total,
           `CAFs` = `Cancer associated fibroblast`/total,
           `Macrophage` = `Macrophage/Monocyte`/total,
           `NK cell` = `NK cell`/total
           ) %>% 
    select("B cell","T cell CD4+","T cell CD8+","CAFs","Macrophage","NK cell")
  
  # Export the results as data.frame 
  assign(paste0(cancer.list[i],".mcp1"), data.frame(mcp.result3))
}


# Save the Rdata 
save.image(file = "celldecon.RData")



