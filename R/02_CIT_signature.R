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
load("data/_raw/CIT_data.Rdata")
load("data/_raw/Bordet.rdata")
GBM_expr <- read_tsv("data/_raw/TCGA.GBM.sampleMAP_HiSeqV2")
GBM_clinical <- read_tsv("data/_raw/TCGA.GBM.sampleMAP_GBM_clinicalMatrix")
GBM_genes <- read_table("data/_raw/gbm_subtype_genes.txt")
source("R/99_func_file.R")

#CIT_full <- probe_to_gene("CIT", "max")
#CIT_full <- CIT_full %>% column_to_rownames(var = "Gene.Symbol")

# Calculate which genes are significantly different for each subtype

save(significant_genes(CIT_full, CIT_classes), file = "data/CIT_signif_subtype_genes.Rdata")
save(significant_genes(GBM_expr, GBM_clinical$GeneExp_Subtype), file = "data/GBM_signif_subtype_genes.Rdata")

save(signatures, file = "data/signatures.Rdata")
load("data/signatures.Rdata")

for (class in unique(CIT_classes)) {
  signa <- data.frame(signatures[[class]])
  colnames(signa) <- c("Probe.Set.ID")
  
  df_joined <- Bordet_annot %>%
    select(Probe.Set.ID, Gene.Symbol) %>% 
    mutate(Gene.Symbol = str_extract(Gene.Symbol, "^\\S+")) %>%
    inner_join(signa, by = "Probe.Set.ID") %>% 
    filter(Gene.Symbol != "---")
  
  file_name <- sprintf("data/signature_CIT_%s_genes.txt", class)
  write.table(df_joined$Gene.Symbol, file = file_name, quote = F, row.names = F, col.names = F)
}

df <- data.frame(performances)
colnames(df) <- c("Subtype", "Sig_len", "Perf")
df$Perf <- as.numeric(as.character(df$Perf))
df$Sig_len <- as.numeric(as.character(df$Sig_len))

df %>% 
  ggplot(aes(x = Sig_len,
             y = Perf,
             color = Subtype)) +
  geom_line() + 
  labs(title = "Performance of variation in subtypes' signature length", x = "Number of genes in signatures", y = "Performance")


# For each sample, the signature with the highest enrichment corresponds to the subtype you assign to the given sample
name_max <- function(column) {
  subtype <- names(column)[which.max(column)]
  return(subtype)
}

enrich <- c()
for (class in c("normL", "lumA", "lumB", "mApo", "basL", "lumC")) {
  signa <- signatures[[class]]
  
  enrich <- rbind(enrich, gsva(CIT_full, 
                               signa,
                               method="ssgsea", 
                               ssgsea.norm = F))
  
}
rownames(enrich) <- c("normL", "lumA", "lumB", "mApo", "basL", "lumC")

pred_class <- apply(enrich, 2, name_max)

mean(CIT_classes == pred_class)

