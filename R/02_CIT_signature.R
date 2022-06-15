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
# load("data/CIT_signif_subtype_genes.Rdata")


## Calculate foldchange signatures
FC_lists <- c()
for (class in unique(CIT_classes)) {
  results <- signif_genes[[class]]
  
  FC_list <- c()
  
  for (probe in names(results)) {
    subtype_median <- median(CIT_full[probe, CIT_classes == class])
    rest_median <- median(CIT_full[probe, CIT_classes != class])
    
    FC_list <- rbind(FC_list, c(probe, log2(subtype_median/rest_median)))
    
  }
  colnames(FC_list) <- c("probe", "FC")
  FC_list <- data.frame(FC_list)
  signa <- arrange(data.frame(FC_list), desc(FC))
  
  file_name <- sprintf("data/ranked_CIT_%s.txt", class)
  
  write(signa$probe, file = file_name)
  
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

for (class in unique(CIT_classes)) {
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

save(signatures, file = "data/signatures.Rdata")

for (class in unique(CIT_classes)) {
  file_name <- sprintf("data/FC_probe_signatures_CIT_%s.txt", class)
  write(signatures[[class]], file = file_name)
  
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
  labs(title = class)


enrich <- c()
for (class in unique(CIT_classes)) {
  signa <- signatures[[class]]
  
  enrich <- rbind(enrich, gsva(CIT_full, 
                               signa,
                               method="ssgsea", 
                               ssgsea.norm = T))
  
}
rownames(enrich) <- unique(CIT_classes)


# For each sample, the signature with the highst enrichment corresponds to the subtype you assign to the given sample
name_max <- function(column) {
  subtype <- names(column)[which.max(column)]
  return(subtype)
}

pred_class <- apply(enrich, 2, name_max)

mean(CIT_classes == pred_class)

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
