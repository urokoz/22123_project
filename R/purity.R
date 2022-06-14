## use the estimate to calculate the enrichment score a stromal 
## content gene set and an immune infiltrate gene set in agiven sample

# download the package 

library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
library('estimate')
library('tidyverse')

prob_extimate <- read.table('data/_raw/CIT_table.txt')

# transform the file into a GCT format 
load("data/_raw/CIT_data.Rdata")
df <- rownames_to_column(data.frame(CIT_full),var = 'NAME')
df$Description <- NA
df <- df %>% relocate(Description, .after = NAME)

write("#1.2", file = "data/_raw/CIT_table.gct")
write.table("54675", file = "data/_raw/CIT_table.gct", append = T, eol = "\t", row.names = F, col.names = F, quote = F)
write.table("355", file = "data/_raw/CIT_table.gct", append = T, row.names = F, col.names = F, quote = F)
write.table(df, file = "data/_raw/CIT_table.gct", append = T,sep = '\t', row.names = F, col.names = T, quote = F)

# calculate the score 
estimateScore(input.ds = 'data/_raw/CIT_table.gct',
              output.ds = 'data/_raw/CIT_scores.gct',
              platform = c("affymetrix"))

# read and tidy the output file


