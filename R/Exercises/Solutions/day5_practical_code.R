#### Course 22123: Computational precision medicine
#### Day 5 practical: Predicting response to checkpoint inhibitor therapy
#### By Lars Ronn Olsen
#### 09/06/2022
#### Technical University of Denmark

## Load necessary packages
library(ggplot2)


## Load the data
# expr contains normalized tpm values as described in the paper by Braun et al.
# pheno contains treatment and response data of which you primarily need to consider the columns:
# RNA_ID: sample names. Matches the column names in expr
# ORR: Objective response rate. CRPR: complete response/partial response, PR: partial response, SD: incomplete response/stable disease, PD: progressive disease, NE: inevaluable
load("~/Dropbox/Teaching/courses/Computational precision medine/Exercises/day 5/day5data.Rdata")


## Inspect tumor mutational burden (TMB) in "CRPR", "PR", "SD", "PD"
pheno$ORR <- factor(pheno$ORR, levels=c("CRPR", "PR", "SD", "PD", "NE"))
ggplot(pheno, aes(x = TMB_Counts, y = ORR)) +
  geom_boxplot() +
  coord_flip()

# QUESTION: Does this work?
# not really...


## Inspect PD-L1 expression (gene name = CD274) in "CRPR", "PR", "SD", "PD"
pheno <- pheno[match(colnames(expr), pheno$RNA_ID),]
df <- data.frame(ORR = pheno$ORR, PDL1 = as.vector(as.matrix(expr[rownames(expr) == "CD274",])))
ggplot(df, aes(x = PDL1, y = ORR)) +
  geom_boxplot() +
  coord_flip()

# QUESTION: Does this work?
# not predictive, but there seems to be a correlation

## Calculate and inspect "T-cell inflamed gene expression profile" in "CRPR", "PR", "SD", "PD"
# Don't worry about the "weighted" sum - just use sum or average. If you read papers by people using this GEP, it's all just average.
TCI_genes <- c("CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT")
TCI_GEP <- colMeans(expr[rownames(expr) %in% TCI_genes,])
df <- data.frame(ORR = pheno$ORR, TCIGEP = TCI_GEP)
ggplot(df, aes(x = TCIGEP, y = ORR)) +
  geom_boxplot() +
  coord_flip()

# QUESTION: Does this work?
# not really

# QUESTION: What do you think about summing genes?
# pretty redimentary approach, but hey, if it works. Unfortunately, it doesn't seem to work for RCC


## Try to combine the 3 metrics somehow and see if you can get better performance - a good place to start is to check if they correlate with each other or whether they explain different variance
df <- data.frame(ORR = pheno$ORR, TCIGEP = TCI_GEP, PDL1 = as.vector(as.matrix(expr[rownames(expr) == "CD274",])), TMB = pheno$TMB_Counts)
ggplot(df, aes(x = TCIGEP, y = TMB)) +
  geom_point(aes(color = ORR)) +
  geom_smooth(method="lm")
ggplot(df, aes(x = TCIGEP, y = PDL1)) +
  geom_point(aes(color = ORR)) + 
  geom_smooth(method = "lm")
ggplot(df, aes(x = PDL1, y = TMB)) +
  geom_point(aes(color = ORR)) + 
  geom_smooth(method = "lm")


## Is response confounded by any of the other variables in the phenotype table?
ggplot(pheno, aes(x = Days_from_TumorSample_Collection_and_Start_of_Trial_Therapy, y = ORR)) +
  geom_boxplot() +
  coord_flip()

# If you find something interesting, think about whether this feature could be confounded with something else
# for example "Days_from_TumorSample_Collection_and_Start_of_Trial_Therapy" - looks promising 
# but can you come up with another explanation as to why responders tends to have a higher value here?
