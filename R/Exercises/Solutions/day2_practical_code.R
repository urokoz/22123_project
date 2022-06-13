#### Course 22123: Computational precision medicine
#### Day 2 practical: working with gene expression data in R
#### By Lars Ronn Olsen
#### 02/06/2022
#### Technical University of Denmark

## Load necessary packages
library("ggplot2")
library("DESeq2")

## Load expression data
gbm_expr <- read.table("~/Dropbox/Teaching/courses/Computational precision medine/Exercises/day 2/day2_data/TCGA.GBM.sampleMap_HiSeqV2", header = TRUE, row.names = 1)
gbm_expr <- round((2^gbm_expr)-1, digits = 0)
colnames(gbm_expr) <- gsub(pattern = "\\.", replacement = "-", x = colnames(gbm_expr))

## Check out the distribution of the expression of a gene or two (or all, if you're brave)
gbm_expr_cpm <- apply(gbm_expr, 2, function(x) x/sum(x)*1000000)
df <- data.frame(gene = gbm_expr_cpm[1,])
ggplot(df, aes(x = gene)) +
  geom_density()

## Load phenotype data
gbm_pheno <- read.table("~/Dropbox/Teaching/courses/Computational precision medine/Exercises/day 2/day2_data/TCGA.GBM.sampleMap_GBM_clinicalMatrix", header = TRUE, row.names = 1, fill = TRUE, sep = "\t")
table(substr(x = rownames(gbm_pheno), start = 14, stop = 15))
gbm_pheno$condition <- rep("GBM", nrow(gbm_pheno))
gbm_pheno$condition[grepl(pattern = "11$", perl = TRUE, rownames(gbm_pheno))] <- "normal"
gbm_pheno$condition[grepl(pattern = "02$", perl = TRUE, rownames(gbm_pheno))] <- "recurrent"
gbm_pheno <- gbm_pheno[rownames(gbm_pheno) %in% colnames(gbm_expr),]
gbm_pheno <- gbm_pheno[match(colnames(gbm_expr), rownames(gbm_pheno)),]
table(rownames(gbm_pheno)==colnames(gbm_expr))

## Optional: save Rdata object for later use
save(gbm_expr, gbm_expr_cpm, gbm_pheno, file = "~/Dropbox/Teaching/courses/Computational precision medine/Exercises/day 2/day2_data/gbm_data.Rdata")
load("~/Dropbox/Teaching/courses/Computational precision medine/Exercises/day 2/day2_data/gbm_data.Rdata")

## Calculate PCA
pca <- prcomp(t(gbm_expr_cpm))

## Check out how much variance is explained
summary(pca)

## Plot PCA
df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])
ggplot(df, aes(x = PC1, y = PC2)) +
  geom_point()

## Plot PCA and color by condition
df$condition <- gbm_pheno$condition
ggplot(df, aes(x = PC1, y = PC2, color=condition)) +
  geom_point()

# QUESTION: Does this look as you would expect? Think about what variance is captured on the different PCs. What happens if you plot other PCs?

## Calculate differentially expressed genes between tumor and normal
gbm_expr_tn <- gbm_expr[, gbm_pheno$condition %in% c("GBM", "normal")]
gbm_pheno_tn <- gbm_pheno[gbm_pheno$condition %in% c("GBM", "normal"), ]
dds <- DESeqDataSetFromMatrix(countData = gbm_expr_tn, colData = gbm_pheno_tn[,colnames(gbm_pheno_tn) %in% c("condition", "GeneExp_Subtype")], design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "normal")
dds <- DESeq(dds)
res <- results(dds)
sum(res$padj < 0.05, na.rm=TRUE)
rownames(res)[res$padj<0.05]

# QUESTION: How well do you think this differential expression analysis will work? Why?

## Subset the expression matrix to differentially expressed genes (padj < 0.05)
gbm_expr_deg <- gbm_expr_cpm[rownames(gbm_expr_cpm) %in% rownames(res)[res$padj<0.00005],]

## Recalculate PCA
pca_deg <- prcomp(t(gbm_expr_deg))

## Check out how much variance is explained
summary(pca_deg)

# QUESTION: Is more or less variance explained? Why? Why not? What did you expect?

## Plot PCA and color by condition
df <- data.frame(PC1 = pca_deg$x[,1], PC2 = pca_deg$x[,2])
df$condition <- gbm_pheno$condition
ggplot(df, aes(x = PC1, y = PC2, color=condition)) +
  geom_point()

## Plot a couple of the top differentially expressed genes
res[order(res$padj),]
df <- data.frame(expr = gbm_expr_cpm[rownames(gbm_expr_cpm) == "FTHL3"])
df$condition <- gbm_pheno$condition
ggplot(df, aes(x = expr, color = condition)) +
  geom_boxplot() +
  coord_flip()

plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

# QUESTION: Are you convinced?

## Plot PCA and color by subtype
pca <- prcomp(t(gbm_expr_cpm))
df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], subtype = gbm_pheno$GeneExp_Subtype)
ggplot(df, aes(x = PC1, y = PC2, color = subtype)) +
  geom_point()

# QUESTION: Does this look as you would expect? Do you recall how the subtypes were defined?

## Subset the expression matrix to verhaak genes
sub_genes <- scan("~/Dropbox/Teaching/courses/Computational precision medine/Exercises/day 2/day2_data/gbm_subtype_genes.txt", what="list")
gbm_expr_cpm_sub <- gbm_expr_cpm[rownames(gbm_expr_cpm) %in% sub_genes, ]

## Run PCA
pca_sub <- prcomp(t(log(gbm_expr_cpm_sub+1)), center = TRUE)
df <- data.frame(PC1 = pca_sub$x[,1], PC2 = pca_sub$x[,2], subtype = gbm_pheno$GeneExp_Subtype)
ggplot(df, aes(x = PC1, y = PC2, color = subtype)) +
  geom_point()
