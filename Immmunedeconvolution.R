# Load packages 
library(dplyr)
library(ggplot2)
library(tidyr)
library(immunedeconv)
library(tibble)
library(biomaRt)
library(TCGAbiolinks) 

################################# Prepare documents ##########################################
setwd("/Users/chensisi/Documents/RNAseq/4_Cell_deconvolution/")

cancer.types <- c("Brain Lower Grade Glioma", "Adrenocortical Cancer","Stomach Adenocarcinoma",
                  "Kidney Clear Cell Carcinoma","Cervical & Endocervical Cancer")

cancer.list <- c("LGG","ACC","STAD","KIRC","CESC")


for (i in 1:length(cancer.types)){
  data <- metadata %>% 
    filter(Project == "TCGA", Sample_type == cancer.types[i]) 
  sampleIDs <- rownames(data)
  filename <- paste0(cancer.list[i], "ids.txt")
  write.table(sampleIDs, filename, col.names = FALSE, row.names = FALSE, quote = FALSE)
}


### Import data 
ensembleID <- read.table("Ensembl.IDs.txt", header = TRUE)
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl", 
                      host = "https://www.ensembl.org/")

files <- c("LGG.expression.txt","ACC.expression.txt","STAD.expression.txt",
           "KIRC.expression.txt","CESC.expression.txt")

for (i in 1:length(files)){
  expression <- read.table(files[i], header = TRUE)
  expression <- cbind(ensembleID, expression)
  
  # Convert ensemblIDs into readable form 
  split <- strsplit(expression$sample, split = "\\.")
  ensemblIDs <- data.frame()
  for (j in 1:nrow(expression)){
    ID <- split[[j]][1]
    ensemblIDs <- rbind(ensemblIDs, ID)
  }
  row.names(expression) <- as.character(ensemblIDs[,1])
  
  # Gene Annotation 
  geneInfo <- getBM(attributes=c("ensembl_gene_id","external_gene_name"),values = row.names(expression),
                    mart=ensembl)
  
  rownames(geneInfo) <- geneInfo$ensembl_gene_id
  
  annotated <- merge(geneInfo, expression, by = 0 )
  annotated2 <- annotated[which(!duplicated(annotated$external_gene_name)),]
  rownames(annotated2) <- annotated2$external_gene_name
  annotated2 <- annotated2[, c(-1,-2,-3,-4)]
  
  # Delog transformation 
  delog <- function(x) {2^x} 
  annotated3 <- as.data.frame(lapply(annotated2, FUN = function(x){sapply(x, FUN = delog)}))
  rownames(annotated3) <- rownames(annotated2)
  
  # File export 
  filename <- paste0(cancer.list[i],"expression.tpm.txt")
  write.table(annotated3, filename, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
}


########################### Cell deconvolution analysis ##########################################
# immunedeconv::deconvolute(gene_expression_matrix, method) 
# methods: quantiseq, timer, cibersort, cibersort_abs, mcp_counter, xcell, epic 

# Test three different cell deconvolution methods using KICH dataset (first 15 samples)\
KICH <- read.table("KICHexpression.tpm.txt", header = TRUE, row.names = 1, sep = "\t")
KICH2 <- as.matrix(KICH[,1:15])
knitr::kable(KICH2[1:5, 1:5])

# EPIC
KICH.epic <- immunedeconv::deconvolute(KICH2, "epic")

# Cibersort
  # set_cibersort_binary("/Users/chensisi/Documents/RNAseq/4_Cell_deconvolution/Cibersort2.R")
  # set_cibersort_mat("/Users/chensisi/Documents/RNAseq/4_Cell_deconvolution/LM22.txt") # PanelA_signature_matrix.txt
  # KICH.cibersort <- immunedeconv::deconvolute(KICH2, "cibersort",perm = 0, absolute = FALSE,arrays = FALSE)

?deconvolute_cibersort(KICH2, arrays = FALSE, absolute = TRUE, abs_method = "sig.score")
source('Cibersort2.R')
KICH.ciber <- CIBERSORT('PanelA_signature_matrix.txt',"KICHexpression.tpm.txt", perm = 100, QN = T) #perm置换次数=1000，QN分位数归一化=TRUE
dim(KICH.ciber)
colnames(KICH.ciber)
KICH.ciber2 <- as.data.frame(KICH.ciber[1:15,])
colnames(KICH.ciber2)
KICH.ciber2 <- KICH.ciber2 %>%
  mutate(`Endothelial cells` = `mv Endothelial cells` + `Endothelial cells`, 
         Macrophages = `Macrophages M1` + `Macrophages M2`) %>%
  dplyr::select("naive B-cells","Fibroblasts","CD4+ T-cells", "CD8+ T-cells","Macrophages",
                "Monocytes","NK cells","Endothelial cells")
  
KICH.ciber3 <- as.data.frame(t(KICH.ciber2))
KICH.ciber3 <- KICH.ciber3 %>% 
  rownames_to_column("cell_type")

# xCell 
  # KICH.xcell <- immunedeconv::deconvolute(KICH2, "xcell")
cell_types = c("B cell", "Cancer associated fibroblast", "T cell CD4+", "T cell CD8+", "Endothelial cell", "Macrophage/Monocyte" , "NK cell")
KICH.xcell1 <- immunedeconv::deconvolute_xcell(KICH2, array = FALSE, expected_cell_types = cell_types)
KICH.xcell1.1 <- as.data.frame(t(KICH.xcell1)) %>% 
  mutate(`CD4+ T-cells` = `CD4+ T-cells` + `Tregs`) %>% 
  dplyr::select("B-cells", "Fibroblasts","CD4+ T-cells","CD8+ T-cells","Endothelial cells","Macrophages","Monocytes","NK cells")

KICH.xcell2 <- as.data.frame(t(KICH.xcell1.1)) %>% 
  rownames_to_column("cell_type")

# quanTISeq 
KICH.quantiseq <- immunedeconv::deconvolute(KICH2, "quantiseq", arrays = FALSE)
KICH.quantiseq2 <- as.data.frame(t(KICH.quantiseq))
colnames(KICH.quantiseq2) <- KICH.quantiseq2[1,]
KICH.quantiseq2 <- KICH.quantiseq2[-1,]
colnames(KICH.quantiseq2)

KICH.quantiseq2 <- KICH.quantiseq2 %>%
  mutate(Macrophages = as.numeric(`Macrophage M1`) + as.numeric(`Macrophage M2`), `CD4+ T cells` = as.numeric(`T cell CD4+ (non-regulatory)`) + as.numeric(`T cell regulatory (Tregs)`)) %>%
  dplyr::select("B cell", "CD4+ T cells","T cell CD8+","Macrophages","uncharacterized cell")

class(KICH.quantiseq2$`Macrophage M1`)

# TIMER  #有bug
KICH.timer <- immunedeconv::deconvolute(KICH2,"timer", indications = "kich")

# MCPcounter 
KICH.mcp <- immunedeconv::deconvolute(KICH2, "mcp_counter")

#### Correlation test between different deconvolution methods ####
# Correlation test between cell deconvolution results
dim(KICH.ciber3)
dim(KICH.epic)
dim(KICH.xcell2)

# Change cell type names to make them consistent 
KICH.ciber3$cell_type
KICH.ciber3$cell_type <- c("B-cells","CAFs", "CD4+ T-cells","CD8+ T-cells","Macrophages","Monocytes","NK cells","Endothelial cells")
KICH.epic$cell_type 
KICH.epic$cell_type <- c("B-cells","CAFs","CD4+ T-cells","CD8+ T-cells","Endothelial cells","Macrophages","NK cells","uncharacterized cells")
KICH.xcell2$cell_type
KICH.xcell2$cell_type <- c("B-cells","CAFs","CD4+ T-cells","CD8+ T-cells","Endothelial cells","Macrophages","Monocytes","NK cells")



#### Data visualization ####
library(dutchmasters)

# Combine all data frames
KICH.ciber4 <- KICH.ciber3
KICH.ciber4$deconvolution <- "CIBERSORT"
KICH.ciber4 <- KICH.ciber4 %>%
  gather(sample, score, -cell_type, -deconvolution)

KICH.epic2 <- KICH.epic
KICH.epic2$deconvolution <- "EPIC"
KICH.epic2 <- KICH.epic %>%
  gather(sample, score, -cell_type, -deconvolution)

KICH.xcell3 <- KICH.xcell2 
KICH.xcell3$deconvolution <- "xCell"
KICH.xcell3 <- KICH.xcell3 %>%
  gather(sample, score, -cell_type, -deconvolution)

KICH.combine <- rbind(KICH.ciber4, KICH.epic2, KICH.xcell3)
unique(KICH.combine$cell_type)
KICH.combine %>%
  filter(cell_type == "B-cells") %>%
  ggplot(aes(x= sample, y = score, color = deconvolution)) + 
  geom_point(size = 2) + 
  facet_wrap(~ cell_type, scales = "free_x", ncol=3)  + 
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1)) 

KICH.combine %>%
  ggplot(aes(x=sample, y=score, color=deconvolution)) +
  geom_point(size=4) +
  facet_wrap(~cell_type, scales="free_x", ncol=3) +
  scale_color_brewer(palette="Paired", guide=FALSE) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

####　Correlation test ####
  #　Data wrangling: transpose the data and set cell_types as column names 
header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

KICH.ciber4 <- as.data.frame(t(KICH.ciber3)) 
KICH.ciber4 <- header.true(KICH.ciber4)
KICH.ciber4$deconvolution <- "CIBERSORT"

KICH.epic4 <- as.data.frame(t(KICH.epic))
KICH.epic4 <- header.true(KICH.epic4)
KICH.epic4$deconvolution <- "EPIC"

KICH.xcell4 <- as.data.frame(t(KICH.xcell2))
KICH.xcell4 <- header.true(KICH.xcell4)
KICH.xcell4$deconvolution <- "xCell"

KICH.xcell5 <- KICH.xcell4 %>% 
  dplyr::select("B-cells","CAFs","CD4+ T-cells","CD8+ T-cells","Macrophages","NK cells","Endothelial cells","deconvolution")
KICH.ciber5 <- KICH.ciber4 %>% 
  dplyr::select("B-cells","CAFs","CD4+ T-cells","CD8+ T-cells","Macrophages","NK cells","Endothelial cells","deconvolution")
KICH.epic5 <- KICH.epic4 %>% 
  dplyr::select("B-cells","CAFs","CD4+ T-cells","CD8+ T-cells","Macrophages","NK cells","Endothelial cells","deconvolution")

cell_types <- colnames(KICH.xcell5)
results <- data.frame()


for (i in 1:8){
  a <- cor.test(as.numeric(KICH.epic4[,i]), as.numeric(KICH.ciber4[,i]), method = "spearman", continuity = TRUE)
  b <- cor.test(as.numeric(KICH.epic4[,i]), as.numeric(KICH.xcell4[,i]), method = "spearman", continuity = TRUE)
  c <- cor.test(as.numeric(KICH.ciber4[,i]), as.numeric(KICH.xcell4[,i]), method = "spearman", continuity = TRUE)
  
  PVAL <- c(a[[3]], b[[3]], c[[3]]) 
  Rho <- c(a[[4]], b[[4]],c[[4]])
  result <- data.frame(PVAL, Rho)
  
  results <- rbind(results, result)
}

cor <- c("EPIC:CIBERSORT","EPIC:xCell","CIBERSORT:xCell") 
correlations <- rep(comparison, 7) 
celltypes  <- rep(colnames(KICH.ciber4)[1:7], each = 3)
cor.results <- data.frame(correlations, celltypes)
cor.results$PVAL <- results$PVAL
cor.results$Rho <- results$Rho




##### Data visualization #####

KICH.epic %>%
  gather(sample, cell_fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=cell_fraction, fill=cell_type)) +
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_dutchmasters(palette = "staalmeesters") +
  scale_x_discrete(limits = rev(levels(KICH.epic))) 

KICH.ciber3 %>%
  gather(sample, cell_fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y= cell_fraction, fill=cell_type)) +
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_dutchmasters(palette = "staalmeesters") +
  scale_x_discrete(limits = rev(levels(KICH.epic)))

KICH.xcell2 %>%
  gather(sample, xCell_enrichment_score, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=xCell_enrichment_score, fill=cell_type)) +  
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_dutchmasters(palette = "staalmeesters") + 
  scale_x_discrete(limits = rev(levels(KICH.xcell2)))

KICH.quantiseq %>%
  gather(sample, cell_fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y= cell_fraction, fill=cell_type)) +
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_dutchmasters(palette = "staalmeesters") +  #scale_fill_brewer(palette="Paired")
  scale_x_discrete(limits = rev(levels(KICH.quantiseq)))

KICH.mcp %>%
  gather(sample, score, -cell_type) %>%
  ggplot(aes(x=sample, y=score, color=cell_type)) +
  geom_point(size=4) +
  facet_wrap(~cell_type, scales="free_x", ncol=3) +
  scale_color_brewer(palette="Paired", guide=FALSE) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))








############# Expected cell types for xCell #########################
# xCell performance can be improved by reducing the
# number of signatures considered. (https://github.com/grst/immunedeconv/issues/1).
# Arguably, many of the xCell signatures do not make sense in
# the context of this benchmark, therefore we limit the signatures
# to the cell types available in the immune datasets.

tmp_quantiseq = KICH.quantiseq %>% map_result_to_celltypes(cell_types, "quantiseq") %>%
  rownames_to_column("cell_type") %>%
  gather("sample", "estimate", -cell_type) %>%
  mutate(method="quanTIseq")


EXPECTED_CELL_TYPES_SC = c("B cell", "Cancer associated fibroblast", "Dendritic cell",
                           "Endothelial cell", "Macrophage/Monocyte", "NK cell",
                           "T cell CD4+ (non-regulatory)", "T cell CD8+",
                           "T cell regulatory (Tregs)")

EXPECTED_CELL_TYPES_FACS = c("B cell", "Dendritic cell", "Monocyte", "NK cell",
                             "T cell CD4+", "T cell CD8+", "T cell")


dataset_racle$expr_mat
KICH2


knitr::kable(dataset_racle$expr_mat[1:5, ])


deconvolution_methods

res_quantiseq = deconvolute(dataset_racle$expr_mat, "quantiseq", tumor = TRUE)


res_quantiseq %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_brewer(palette="Paired") +
  scale_x_discrete(limits = rev(levels(res_quantiseq)))

res_mcp_counter = deconvolute(dataset_racle$expr_mat, "mcp_counter")

res_mcp_counter %>%
  gather(sample, score, -cell_type) %>%
  ggplot(aes(x=sample, y=score, color=cell_type)) +
  geom_point(size=4) +
  facet_wrap(~cell_type, scales="free_x", ncol=3) +
  scale_color_brewer(palette="Paired", guide="none") +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

cell_types = c("B cell", "T cell CD4+", "T cell CD8+", "NK cell")

tmp_quantiseq = res_quantiseq %>% 
  map_result_to_celltypes(cell_types, "quantiseq") %>%
  rownames_to_column("cell_type") %>%
  gather("sample", "estimate", -cell_type) %>%
  mutate(method="quanTIseq")

tmp_mcp_counter = res_mcp_counter %>% map_result_to_celltypes(cell_types, "mcp_counter") %>%
  rownames_to_column("cell_type") %>%
  gather("sample", "estimate", -cell_type) %>%
  mutate(method="MCP-counter")

result = bind_rows(tmp_quantiseq, tmp_mcp_counter) %>%
  inner_join(dataset_racle$ref)


result %>%
  ggplot(aes(x=true_fraction, y=estimate)) +
  geom_point(aes(shape=cell_type, color=cell_type)) +
  facet_wrap(cell_type~method, scales="free_y", ncol = 2) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw()











KICH <- read.table("TCGA-KICH.htseq_fpkm.tsv", header = TRUE, row.names = 1)
LUSC <- read.table("TCGA.LUSC.sampleMap_HiSeqV2", header = TRUE,  row.names = 1)

library(biomaRt)
library(TCGAbiolinks) 

ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl", 
                      mirror = "useast")

geneInfo <- getBM(attributes=c("ensembl_transcript_id","external_gene_name"),values = row.names(KICH),
                  mart=ensembl)

geneInfo <- geneInfo[which(!geneInfo$external_gene_name==""),]
geneInfo <- geneInfo[which(!duplicated(geneInfo$external_gene_name)),]











xcell.KICH %>% 
  gather(samples, cell_fraction, c("TCGA.KN.8421.01","TCGA.KO.8417.01","TCGA.KM.8438.01",
                              "TCGA.KL.8340.01","TCGA.KO.8406.01", "TCGA.KO.8408.01",
                              "TCGA.KM.8443.01", "TCGA.KM.8442.01", "TCGA.KL.8332.11",
                              "TCGA.KL.8327.01")) %>%
  ggplot(aes(x = samples, y = cell_fraction , fill = cell_type)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_y_continuous(labels = scales::percent_format()) + 
  labs(y = "Cell fraction", x = "KICH samples") + 
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1), 
        plot.title = element_text(hjust = 0.5))

xcell2.KICH <- xcell2.KICH %>% 
  select(1:21)

xcell2.KICH %>% 
  gather(samples, cell_fraction, 2:21) %>%
  ggplot(aes(x = samples, y = cell_fraction, fill = Cell_type)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_y_continuous(labels = scales::percent_format()) + 
  labs(y = "Cell fraction", x = "KICH samples") + 
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1), 
        plot.title = element_text(hjust = 0.5))

library(RColorBrewer)
n <- 67
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


# scale_fill_dutchmasters(palette = "staalmeesters")
# scale_fill_igv() or  scale_fill_ucscgb()

devtools::install_github("jrnold/ggthemes")
library(ggthemes)
?ggthemes

dim(xcell.KICH)
display.brewer.all()
devtools::install_github("road2stat/ggsci")
install.packages("ggsci")
library(ggsci)


setwd("/Users/chensisi/Documents/RNAseq/2_Survival analysis /")

write.table(multivariate.results, "multivariate.results_by.ratio.csv", sep =",", col.names = TRUE, row.names = FALSE)
