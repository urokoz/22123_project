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


# Calculate which genes are significantly different for each subtype
signif_genes <- c()

for (class in unique(CIT_classes)) {
  
  results <- c()
  for (probe in rownames(CIT_full)) {
    is_class <- CIT_full[probe, CIT_classes == class]
    rest <- CIT_full[probe, CIT_classes != class]
    
    res <- wilcox.test(x = is_class,
                       y = rest)
    results[[probe]] <- res$p.value
    
  }
  results <- p.adjust(results)
  results <- results[results < 0.05]
  
  signif_genes[[class]] <- results
  
}

save(signif_genes, file = "data/CIT_signif_subtype_genes.Rdata")


## Calculate foldchange signatures
FC_lists <- c()
for (class in unique(CIT_classes)) {
  results <- signif_genes[[class]]
  
  FC_list <- c()
  
  for (probe in names(results)) {
    subtype_median <- median(CIT_full[probe, CIT_classes == class])
    rest_median <- median(CIT_full[probe, CIT_classes != class])
    
    FC_list[[probe]] <- log2(subtype_median/rest_median)
    
  }
  
  signa <- sort(FC_list, decreasing = T)
  
  file_name <- sprintf("data/ranked_CIT_%s.txt", class)
  
  write(names(signa), file = file_name)
  
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
perf_collect <- c()

for (class in unique(CIT_classes)) {
  sig_probes <- FC_lists[[class]]
  
  best_perf <- 0
  conseq_worse <- 0
  best_size <- NULL
  performances <- c()
  for (i in 1:length(sig_probes)) {
    signa <- list(names(sig_probes[1:i]))
    
    perf <- pred_perf(signa, class)
    print(sprintf("%s, %d:  %f", class, i, perf))
    performances[i] <- perf
    
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
  signatures[[class]] <- list(names(sig_probes[1:best_size]))
  perf_collect[[class]] <- performances
}


for (class in unique(CIT_classes)) {
  df <- perf_collect[[class]]
  ggplot(aes(x = length(df), y = df)) +
    geom_line()
  
}



# for each class
#   for i in sig probe 
#     define signature
#     
#     run performace (ssgsea)
#     
#     save performance 
#     
#     save i if best performance
#     
#     if worse 5 x in row 
#       break
#   
#   save performace history
#   save best signature for class
#   
# plot performances 
# 
# run collective performance
