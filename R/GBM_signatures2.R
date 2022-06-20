#### in this section I will calculate and rank the enrichment of different genes in
#### sample level and see if there are some genes are constantly on the top 25% of the
#### specific subtype and in the button of the rest.

## load the expression data
## the CIT_full contains all the expression data for all probes 
## the CIT_annot contains annotions of the 375 probes used for classification  
## the CIT class contains the subtypes of all the samples

mean_expression(GBM_expr, GBM_classes)
## find the overlap of the genes that exist in both the top of the specifical class
## and the buttom of the rest.

Overlap_genes <- function(subtypes, percentage){
  temp_class <- read.table(sprintf('data/top_or_buttom_25/%s_class_GBM.txt',subtypes))
  temp_rest  <- read.table(sprintf('data/top_or_buttom_25/%s_rest_GBM.txt',subtypes))
  
  index_25 <- as.integer(nrow(temp_class)*percentage)
  index_75 <- as.integer(nrow(temp_rest)*(1-percentage))
  top_class_list <- rownames(temp_class)[1:index_25]
  bot_rest_list <- rownames(temp_rest)[index_75:nrow(temp_rest)]
  overlap <- top_class_list[top_class_list %in% bot_rest_list]
  {write(overlap, file =sprintf('data/top_or_buttom_25/0.35/%s_signatures_GBM.txt',subtypes))}}


percentage <- 0.35
# test
for(subtypes in unique(GBM_classes)){Overlap_genes(subtypes, percentage)}
