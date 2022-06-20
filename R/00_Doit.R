library('tidyverse')
library('broom')
library('GSVA')
library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
library('estimate')
source("R/99_func_file.R")

#pick pipeline
#CIT
load("data/_raw/CIT_data.Rdata")
load("data/_raw/Bordet.rdata")
source("R/02_CIT_signature.R")
source("R/03_CIT_signature1.R")
source("R/04_CIT_signature3.R")
source("R/05_purity.R")


#GBM
GBM_expr <- read_tsv("data/_raw/TCGA.GBM.sampleMAP_HiSeqV2")
GBM_clinical <- read_tsv("data/_raw/TCGA.GBM.sampleMAP_GBM_clinicalMatrix")
GBM_genes <- read_table("data/_raw/gbm_subtype_genes.txt")
source("R/01_load.R")
source("R/GBM_signatures1.R")
source("R/GBM_signatures2.R")
source("R/GBM_signatures3.R")
source("R/purity_GBM.R")
