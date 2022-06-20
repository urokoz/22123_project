#### Course 22123: Computational precision medicine
#### Project work
#### By Mathias Rahbek-Borre, Yu Liu and Joachim Breitenstein
#### 20/06/2022
#### Technical University of Denmark

## Load necessary packages
library('GSVA')
library("tidyverse")

## Load expression data (NOTE: this is microarray data, not RNA-seq!)
# CIT_full contains the expression of all probes for all samples
# CIT_annot contains the annotation for the 375 probes used for classification
# CIT_class contains the subtypes
#load("data/_raw/CIT_data.Rdata")
#load("data/_raw/Bordet.rdata")
#GBM_expr <- read_tsv("data/_raw/TCGA.GBM.sampleMAP_HiSeqV2")
#GBM_clinical <- read_tsv("data/_raw/TCGA.GBM.sampleMAP_GBM_clinicalMatrix")
#GBM_genes <- read_table("data/_raw/gbm_subtype_genes.txt")
#source("R/99_func_file.R")

#CIT_full <- probe_to_gene("CIT", "max")
#CIT_full <- CIT_full %>% column_to_rownames(var = "Gene.Symbol")



# Calculate which genes are significantly different for each subtype
signif_genes_GBM <- significant_genes(GBM_expr, GBM_classes)
save(signif_genes_GBM, file = "data/mann_whitney_u_test/GBM_signif_subtype_genes.Rdata")

FC_GBM <- FC_calc(GBM_expr, signif_genes_GBM, GBM_classes)
save(FC_GBM, file = "data/mann_whitney_u_test/GBM_FC.Rdata")

signatures_GBM <- calc_signatures(GBM_expr, FC_GBM, GBM_classes)

performances <- signatures_GBM$performances
signatures <- signatures_GBM$signatures

save(signatures, file = "data/mann_whitney_u_test/GBM_signatures.Rdata")


for (class in unique(GBM_classes)) {
  signa <- data.frame(signatures[[class]])
  colnames(signa) <- c("Gene.Symbol")

  file_name <- sprintf("data/mann_whitney_u_test/signature_GBM_%s_genes.txt", class)
  write.table(signa, file = file_name, quote = F, row.names = F, col.names = F)
}
# 
# df <- data.frame(performances)
# colnames(df) <- c("Subtype", "Sig_len", "Perf")
# df$Perf <- as.numeric(as.character(df$Perf))
# df$Sig_len <- as.numeric(as.character(df$Sig_len))
# 
# df %>%
#   ggplot(aes(x = Sig_len,
#              y = Perf,
#              color = Subtype)) +
#   geom_line() +
#   labs(title = "Performance of variation in subtypes' signature length", x = "Number of genes in signatures", y = "Performance")


# # For each sample, the signature with the highest enrichment corresponds to the subtype you assign to the given sample
# name_max <- function(column) {
#   subtype <- names(column)[which.max(column)]
#   return(subtype)
# }
# 
# pred_classes <- function(df, classes) {
#   enrich <- c()
#   for (class in unique(classes)) {
#     signa <- signatures[[class]]
#     
#     enrich <- rbind(enrich, gsva(df, 
#                                  signa,
#                                  method="ssgsea", 
#                                  ssgsea.norm = F))
#     
#   }
#   rownames(enrich) <- unique(classes)
#   
#   pred_class <- apply(enrich, 2, name_max)
#   
#   mean(classes == pred_class)
# }
# 
# 
