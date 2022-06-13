#### Course 22123: Computational precision medicine
#### Day 4 practical: Dealing with unwanted variance when subtyping
#### By Lars Ronn Olsen
#### 07/06/2022
#### Technical University of Denmark

## Load necessary packages
library(ggplot2)
library(GSVA)
library("tidyverse")


## Load expression data (NOTE: this is microarray data, not RNA-seq!)
# CIT_full contains the expression of all probes for all samples
# CIT_annot contains the annotation for the 375 probes used for classification
# CIT_class contains the subtypes
load("data/CIT_data.Rdata")


## If you didn't do this yesterday, try to visualize the six subtypes using PCA


pca_fit <- prcomp(data.frame(t(CIT_full)))

df <- data.frame(PC1 = pca_fit$x[,1], PC2 = pca_fit$x[,2])
df$subtype <- CIT_classes
df %>% ggplot(aes(x = PC1,
                  y = PC2,
                  color = subtype)) +
  geom_point()


## Define a signature for each subtype 
# Using CIT_full, run a Mann-Whitney U test for each probe, for samples in each subtype vs the rest of the subtypes combined
# For each comparison, save probes that are significantly differentially expressed (remember multiple testing correction! Use the function "p.adjust")
# For significantly differentially expressed probes, calculate the log2 fold change (log2(median(samples-in-subtype)/median(samples-in-rest-of-the-subtypes)))
# Pick the highest upregulated probes (either by a log2fc threshold like >1, or by some reasonable size - maybe 40 probes?)
# This is your signature for this subtype

signatures <- c()

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
  
  FC_list <- c()
  
  for (probe in names(results)) {
    subtype_median <- median(CIT_full[probe, CIT_classes == class])
    rest_median <- median(CIT_full[probe, CIT_classes != class])
    
    FC_list[[probe]] <- log2(subtype_median/rest_median)
    
  }
  
  signa <- sort(FC_list, decreasing = T)[1:40]
  
  signatures[[class]] <- signa
  
}

## For each sample, use ssGSEA (GSVA) to calculate the enrichment of each signature
# "?GSVA"..

enrich <- c()
for (class in unique(CIT_classes)) {
  signa <- list(names(signatures[[class]]))
  
  enrich <- rbind(enrich, gsva(CIT_full, 
                 signa,
                 method="ssgsea"))
  
}
rownames(enrich) <- unique(CIT_classes)


# For each sample, the signature with the highst enrichment corresponds to the subtype you assign to the given sample
name_max <- function(column) {
  subtype <- names(column)[which.max(column)]
  return(subtype)
}

pred_class <- apply(enrich, 2, name_max)

mean(CIT_classes == pred_class)

## If time: load the paired microarray and RNA-seq data set


## Plot the counts vs probe intensities

# QUESTION: Does it correlate?


## Plot the tpm vs probe intensities

# QUESTION: Does it correlate better?


## Convert the tpm and probe intensities to ranks and plot them

# QUESTION: How about the correlation now?


## Try to classify each sample in the microarray dataset and the RNAseq (TPM) dataset using the signature enrichment approach

# QUESTION: Do you get the same classses?
