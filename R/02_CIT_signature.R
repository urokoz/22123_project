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
source("R/99_func_file.R")

#CIT_full <- probe_to_gene("CIT", "max")
#CIT_full <- CIT_full %>% column_to_rownames(var = "Gene.Symbol")

# Calculate which genes are significantly different for each subtype

signif_genes <- c()

for (class in unique(CIT_classes)) {
  
  results <- c()
  for (probe in rownames(CIT_full)) {
    is_class <- CIT_full[probe, CIT_classes == class]
    rest <- CIT_full[probe, CIT_classes != class]
    res <- wilcox.test(x = as.numeric(is_class),
                       y = as.numeric(rest))
    results[[probe]] <- res$p.value
    
  }
  results <- p.adjust(results)
  results <- results[results < 0.05]
  
  signif_genes[[class]] <- results
  
}

save(signif_genes, file = "data/mann_whitney_u_test/CIT_signif_subtype_genes.Rdata")
# load("data/mann_whitney_u_test/CIT_signif_subtype_genes.Rdata")


## Calculate foldchange signatures
FC_lists <- c()
for (class in unique(CIT_classes)) {
  results <- signif_genes[[class]]
  
  FC_list <- c()
  
  for (probe in names(results)) {
    subtype_median <- median(as.numeric(CIT_full[probe, CIT_classes == class]))
    rest_median <- median(as.numeric(CIT_full[probe, CIT_classes != class]))
    
    FC_list <- rbind(FC_list, c(probe, log2(subtype_median/rest_median)))
    
  }
  colnames(FC_list) <- c("probe", "FC")
  FC_list <- data.frame(FC_list)
  signa <- arrange(data.frame(FC_list), desc(FC))
  
  file_name <- sprintf("data/mann_whitney_u_test/ranked_CIT_%s_probes.txt", class)
  
  write.table(signa$probe, file = file_name, quote = F, row.names = F, col.names = F)
  
  FC_lists[[class]] <- signa
  
}

pred_perf <- function(signa, class) {
  enrich <- gsva(CIT_full,
                 signa,
                 method="ssgsea",
                 ssgsea.norm = F,
                 verbose=F)
  return(ks.test(enrich[CIT_classes == class], enrich[CIT_classes != class])$statistic)
}

signatures <- c()
performances <- c()

for (class in c("lumB", "lumC", "basL")) {
  sig_probes <- FC_lists[[class]]$probe
  
  best_perf <- 0
  conseq_worse <- 0
  best_size <- NULL
  for (i in 1:length(sig_probes)) {
    signa <- list(sig_probes[1:i])
    
    perf <- pred_perf(signa, class)
    print(sprintf("%s, %d:  %f", class, i, perf))
    performances <- rbind(performances, t(c(class, i, perf)))
    
    if (perf <= best_perf) {
      conseq_worse <- conseq_worse + 1
      if (conseq_worse == 25) {
        break
      }
    } else {
      conseq_worse <- 0
      best_perf <- perf
      best_size <- i
    }
  }
  signatures[[class]] <- list(sig_probes[1:best_size])
}

save(signatures, file = "data/mann_whitney_u_test/signatures.Rdata")
load("data/mann_whitney_u_test/signatures.Rdata")

for (class in unique(CIT_classes)) {
  signa <- data.frame(signatures[[class]])
  colnames(signa) <- c("Probe.Set.ID")
  
  df_joined <- Bordet_annot %>%
  select(Probe.Set.ID, Gene.Symbol) %>% 
  mutate(Gene.Symbol = str_extract(Gene.Symbol, "^\\S+")) %>%
  inner_join(signa, by = "Probe.Set.ID") %>% 
  filter(Gene.Symbol != "---")
  
  file_name <- sprintf("data/mann_whitney_u_test/signature_CIT_%s_genes.txt", class)
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


