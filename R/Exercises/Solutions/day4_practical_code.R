#### Course 22123: Computational precision medicine
#### Day 4 practical: Dealing with unwanted variance when subtyping
#### By Lars Ronn Olsen
#### 08/06/2022
#### Technical University of Denmark

## Load necessary packages
library(ggplot2)
library(GSVA)


## Load expression data (NOTE: this is microarray data, not RNA-seq!)
# CIT_full contains the expression of all probes for all samples
# CIT_annot contains the annotation for the 375 probes used for classification
# CIT_class contains the subtypes
load("~/Dropbox/Teaching/courses/Computational precision medine/Exercises/day 3/CIT_data.Rdata")


## If you didn't do this yesterday, try to visualize the six subtypes using PCA
# Run PCA on all samples, all probes
pca <- prcomp(t(CIT_full), scale = TRUE)
df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], subtype = CIT_classes)
ggplot(df, aes(x = PC1, y = PC2, color = subtype)) +
  geom_point()
# Run PCA on all samples, only CIT probes
CIT_sub <- CIT_full[rownames(CIT_full) %in% CIT_annot$Probe.Set.ID,]
pca <- prcomp(t(CIT_sub), scale = TRUE)
df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], subtype = CIT_classes)
ggplot(df, aes(x = PC1, y = PC2, color = subtype)) +
  geom_point()


## Define a signature for each subtype 
# Using CIT_full, run a Mann-Whitney U test for each probe, for samples in each subtype vs the rest of the subtypes combined
# For each comparison, save probes that are significantly differentially expressed (remember multiple testing correction! Use the function "p.adjust")
# For significantly differentially expressed probes, calculate the log2 fold change (log2(median(samples-in-subtype)/median(samples-in-rest-of-the-subtypes)))
# Pick the highest upregulated probes (either by a log2fc threshold like >1, or by some reasonable size - maybe 40 probes?)
# This is your signature for this subtype

# Make an empty list to store your signatures in
probe_signatures <- list()
# Make a function to do a MWW U test
row_mww <- function(x) {
  subtype <- x[CIT_classes == i]
  rest <- x[!CIT_classes == i]
  res <- wilcox.test(subtype, rest)
  return(res$p.value)
}
# Make a function to calculate log2 fc
row_fc <- function(x) {
  subtype <- x[CIT_classes == i]
  rest <- x[!CIT_classes == i]
  res <- log2(median(subtype)/median(rest))
  return(res)
}
# Run the comparisons for each CIT class
for (i in unique(CIT_classes)) {
  pvals <- apply(CIT_full, 1, FUN = row_mww)
  padj <- p.adjust(pvals, method = "BY")
  sig_probes <- rownames(CIT_full)[padj<0.05]
  subtype_sig <- CIT_full[rownames(CIT_full) %in% sig_probes,]
  log2fc <- apply(subtype_sig, 1, FUN = row_fc)
  probe_signatures[[i]] <- names(log2fc[order(log2fc, decreasing = TRUE)][1:50])
  print(i)
}


## Bonus figure: are we capturing our subtypes?
library(ComplexHeatmap)
column_ha = HeatmapAnnotation(subtype = CIT_classes)
Heatmap(scale(CIT_full[rownames(CIT_full) %in% unname(unlist(probe_signatures)),]), show_column_names = FALSE, show_row_names = FALSE, top_annotation = column_ha)

# QUESTION: Does the subtypes cluster together? If no, what do you think went wrong?


## For each sample, use ssGSEA (GSVA) to calculate the enrichment of each signature
# "?gsva"
# For each sample, the signature with the highest enrichment corresponds to the subtype you assign to the given sample
scores <- gsva(CIT_full, probe_signatures, method = "ssgsea", ssgsea.norm = FALSE)
table(CIT_classes == rownames(scores)[apply(scores, 2, which.max)])


## If time: load the paired microarray and RNA-seq data set
load("~/Dropbox/Teaching/courses/Computational precision medine/Exercises/day 4/Bordet.rdata")


## Plot the counts vs probe intensities (for a sample)
df <- data.frame(microarray = Bordet_array[,1], RNAseq = log(Bordet_RNA_count[,1]+1))
ggplot(df, aes(x = microarray, y = RNAseq)) +
  geom_point() +
  geom_smooth()

# QUESTION: Does it correlate?
cor(x = Bordet_array[,1], y = log(Bordet_RNA_count[,1]+1), method = "pearson")
# Yes, kind of


## Plot the tpm vs probe intensities
df <- data.frame(microarray = Bordet_array[,1], RNAseq = log(Bordet_RNA_tpm[,1]+1))
ggplot(df, aes(x = microarray, y = RNAseq)) +
  geom_point() +
  geom_smooth()

# QUESTION: Does it correlate better?
cor(x = Bordet_array[,1], y = log(Bordet_RNA_tpm[,1]+1), method = "pearson")
# Similar


## Convert the tpm and probe intensities to ranks and plot them
df <- data.frame(microarray = rank(Bordet_array[,1]), RNAseq = rank(Bordet_RNA_tpm[,1]))
ggplot(df, aes(x = microarray, y = RNAseq)) +
  geom_point() +
  geom_smooth()

# QUESTION: How about the correlation now?
cor(x = rank(Bordet_array[,1]), y = rank(Bordet_RNA_tpm[,1]), method = "pearson")
# Still good!


## Try to classify each sample in the microarray dataset and the RNAseq (TPM) dataset using the signature enrichment approach
microarray_scores <- gsva(Bordet_array, probe_signatures, method = "ssgsea", ssgsea.norm = FALSE)
microarray_subtypes <- rownames(microarray_scores)[apply(scores, 2, which.max)]
RNAseq_scores <- gsva(Bordet_array, probe_signatures, method = "ssgsea", ssgsea.norm = FALSE)
RNAseq_subtypes <- rownames(RNAseq_scores)[apply(scores, 2, which.max)]

# QUESTION: Do you get the same classses?
table(microarray_subtypes == RNAseq_subtypes)
# Yay!


## Try to classify each samples microarray dataset and the RNAseq (TPM) dataset using the distance to centroid approach
## Note: you will need to convert all vectors to rank - including the training set and your centroids, and use the Kendall Tau distance metric to centroids