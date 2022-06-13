#### Course 22123: Computational precision medicine
#### Day 5 practical: Predicting response to checkpoint inhibitor therapy
#### By Lars Ronn Olsen
#### 09/06/2022
#### Technical University of Denmark

## Load necessary packages
library(ggplot2)
library("tidyverse")
library(ggpubr)

## Load the data
# expr contains normalized tpm values as described in the paper by Braun et al.
# pheno contains treatment and response data of which you primarily need to consider the columns:
# RNA_ID: sample names. Matches the column names in expr
# ORR: Objective response rate. CRPR: complete response/partial response, PR: partial response, SD: incomplete response/stable disease, PD: progressive disease, NE: inevaluable
load("data/day5data.Rdata")

## Inspect tumor mutational burden (TMB) in "CRPR", "PR", "SD", "PD"
pheno %>% 
  filter(ORR != "NE") %>% 
  ggplot(aes(x = ORR,
             y = TMB_Counts,
             fill = ORR)) +
  geom_violin() + 
  geom_boxplot(width = 0.05,
               color = "black",
               fill = "white",
               outlier.shape = NA)

pheno %>% 
  filter(ORR != "NE") %>% 
  count(ORR)

# QUESTION: Does this work?
# Not really. while there are some variation between the groups, they do not seem particularly different. Especially considering that the groups are quite small.

## Inspect PD-L1 expression (gene name = CD274) in "CRPR", "PR", "SD", "PD"
ORR_df <- data.frame(pheno[, c("RNA_ID", "ORR")]) 
CD274_df <- pivot_longer(expr["CD274", ],everything(), names_to = "RNA_ID", values_to = "CD274")

df <- inner_join(ORR_df, CD274_df, by = "RNA_ID")

df %>% 
  filter(ORR != "NE") %>%  
  ggplot(aes(x = ORR,
             y = CD274,
             fill = ORR)) +
  geom_violin() + 
  geom_boxplot(width = 0.05,
               color = "black",
               fill = "white",
               outlier.shape = NA)
  
# QUESTION: Does this work?
# It doesn't really separate the ORR groupings. The expressions are quite similar for the cohorts.

## Calculate and inspect "T-cell inflamed gene expression profile" in "CRPR", "PR", "SD", "PD"
# Don't worry about the "weighted" sum - just use sum or average. If you read papers by people using this GEP, it's all just average.
GEP_genes <- c("CCL5","CD27","CD274","CD276","CD8A","CMKLR1","CXCL9","CXCR6","HLA-DQA1","HLA-DRB1","HLA-E","IDO1","LAG3","NKG7","PDCD1LG2","PSMB10","STAT1","TIGIT")

GEP_df <- expr[GEP_genes, ] %>% 
  drop_na() %>% 
  rownames_to_column(var = "Genes") %>% 
  pivot_longer(!Genes, 
               names_to = "RNA_ID", 
               values_to = "Expression")

GEP_df <- inner_join(ORR_df, GEP_df, by = "RNA_ID") %>% 
  filter(ORR != "NE")

GEP_df %>% 
  group_by(ORR, RNA_ID) %>% 
  summarise(mean = mean(Expression)) %>%  
  ggplot(aes(x = ORR,
             y = mean,
             fill = ORR)) +
  geom_violin() + 
  geom_boxplot(width = 0.05,
               color = "black",
               fill = "white",
               outlier.shape = NA)




# QUESTION: Does this work?
# Not really. The plots look like the last.

# QUESTION: What do you think about summing genes?
# Generally it seems like a silly idea. It makes it very difficult to see how the genes in that signature is affected. 
# Maybe some are upregulated and others are turned off. And this could differ from patient to patient.


## Try to combine the 3 metrics somehow and see if you can get better performance - a good place to start is to check if they correlate with each other or whether they explain different variance


## Is response confounded by any of the other variables in the phenotype table?
pheno %>% 
  filter(ORR != "NE") %>% 
  ggplot(aes(x = ORR,
             y = Age,
             fill = ORR)) +
  geom_violin() + 
  geom_boxplot(width = 0.05,
               color = "black",
               fill = "white",
               outlier.shape = NA)

pheno %>% 
  filter(ORR != "NE") %>% 
  ggplot(aes(x = ORR,
             y = Days_from_TumorSample_Collection_and_Start_of_Trial_Therapy,
             fill = ORR)) +
  geom_violin() + 
  geom_boxplot(width = 0.05,
               color = "black",
               fill = "white",
               outlier.shape = NA)
