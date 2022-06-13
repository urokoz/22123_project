#### Course 22123: Computational precision medicine
#### Day 2 practical: working with gene expression data in R
#### By Mathias Rahbek-Borre
#### 02/06/2022
#### Technical University of Denmark

## Load necessary packages
library("tidyverse")
library("broom")


## Load expression data
gbm_expr <- read_tsv("data/TCGA.GBM.sampleMap_HiSeqV2.gz") %>% 
  pivot_longer(!sample, 
               names_to = "sampleID",
               values_to = "counts") %>%
  mutate(counts = round((2^counts)-1)) %>% 
  mutate(counts = (counts / sum(counts))*1000000) %>% 
  pivot_wider(names_from = sample, values_from = counts)
  

## Check out the distribution of the expression of a gene or two (or all, if you're brave)
gbm_expr %>% 
  ggplot(aes(x = RTN4RL2)) +
  geom_density()

## Load phenotype data
gbm_pheno <- read_tsv("data/TCGA.GBM.sampleMap_GBM_clinicalMatrix")

gbm_all <- inner_join(gbm_expr, gbm_pheno, by = "sampleID")
## Optional: save Rdata object for later use


## Calculate PCA
pca_fit <- gbm_expr %>% 
  select(!sampleID) %>% 
  prcomp()

## Check out how much variance is explained
pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  filter(percent > 0.025) %>%
  ggplot(mapping = aes(x = PC,
                       y = percent)) +
  geom_col(fill = "#56B4E9", 
           alpha = 0.8) +
  scale_x_continuous(breaks = 1:5) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.01))) +
  theme_minimal(base_size = 19) +
  theme(plot.background = element_rect(colour = "black", 
                                       fill = NA,
                                       size = 1)) +
  labs(y = "Percent", 
       title = "Variance explained by each principle component")

## Plot PCA
pca_fit %>%
  augment(gbm_all) %>% # add original dataset back in
  ggplot(mapping = aes(x = .fittedPC1, 
                       y = .fittedPC2)) +
  geom_point(size = 1.5)

## Plot PCA and color by condition
pca_fit %>%
  augment(gbm_all) %>% # add original dataset back in
  ggplot(mapping = aes(x = .fittedPC1, 
                       y = .fittedPC2, 
                       color = GeneExp_Subtype)) +
  geom_point(size = 1.5)

## Calculate differentially expressed genes between tumor and normal


## Subset the expression matrix to differentially expressed genes (padj < 0.05)


## Recalculate PCA


## Check out how much variance is explained


## Plot PCA and color by condition


## Plot a couple of the top differentially expressed genes


## Plot PCA and color by subtype


## Subset the expression matrix to verhaak genes


## Run PCA

