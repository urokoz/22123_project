## use the estimate to calculate the enrichment score a stromal 
## content gene set and an immune infiltrate gene set in agiven sample

# download the package 

library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
library('estimate')
library('tidyverse')
source("R/99_func_file.R")

# transform the file into a GCT format 
load("data/_raw/CIT_data.Rdata")

df_master <- probe_to_gene("CIT", "max")
df <- df_master[2:nrow(df_master),]

names(df)[names(df) == "Gene.Symbol"] <- "NAME"

df <- column_to_rownames(df, var = "NAME")

write.table(df, file = "data/_raw/CIT_samples.txt", sep = "\t", quote = F)

filterCommonGenes(input.f="data/_raw/CIT_samples.txt", output.f='data/_raw/CIT_genes.gct', id="GeneSymbol")

# calculate the score 
estimateScore(input.ds = 'data/_raw/CIT_genes.gct',
              output.ds = 'data/_raw/CIT_scores.gct',
              platform = c("affymetrix"))


write.table(df$NAME, file = "regex.test", row.names = F, col.names = F, quote = F)
# read and tidy the output file

