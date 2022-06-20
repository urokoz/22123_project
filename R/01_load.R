#### Course 22123: Computational precision medicine
#### Project work
#### By Mathias Rahbek-Borre, Yu Liu and Joachim Breitenstein
#### 20/06/2022
#### Technical University of Denmark

## Load necessary packages
library(GSVA)
library("tidyverse")

## Load expression data (NOTE: this is microarray data, not RNA-seq!)
# CIT_full contains the expression of all probes for all samples
# CIT_annot contains the annotation for the 375 probes used for classification
# CIT_class contains the subtypes
load("data/_raw/CIT_data.Rdata")
load("data/_raw/Bordet.rdata")
GBM_expr <- read_tsv("data/_raw/TCGA.GBM.sampleMAP_HiSeqV2")
GBM_clinical <- read_tsv("data/_raw/TCGA.GBM.sampleMAP_GBM_clinicalMatrix")
GBM_genes <- read_table("data/_raw/gbm_subtype_genes.txt")
source("R/99_func_file.R")



GBM_transposed <- GBM_expr %>% 
  as.data.frame() %>% 
  t()

colnames(GBM_transposed) <- GBM_transposed[1, ]
GBM_transposed <- GBM_transposed[-1, ] %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sampleID")

GBM_classes <- GBM_transposed %>% 
  inner_join(GBM_clinical, by = "sampleID") %>% 
  select(GeneExp_Subtype) 

GBM_classes <- GBM_classes$GeneExp_Subtype

GBM_expr <- column_to_rownames(GBM_expr, var = "sample")

GBM_classes[is.na(GBM_classes)] <- "Unknown"


#test_classes <- data.frame(samples = GBM_clinical$sampleID, subtypes = GBM_clinical$GeneExp_Subtype)

#test_classes <- column_to_rownames(test_classes, var = "samples")

#test_classes2 <- as.character(test_classes[colnames(GBM_expr),])
