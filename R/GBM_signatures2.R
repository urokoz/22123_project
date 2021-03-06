#### in this section I will calculate and rank the enrichment of different genes in
#### sample level and see if there are some genes are constantly on the top 25% of the
#### specific subtype and in the button of the rest.

## load the expression data
## the CIT_full contains all the expression data for all probes 
## the CIT_annot contains annotions of the 375 probes used for classification  
## the CIT class contains the subtypes of all the samples

mean_expression <- function(df, classes){
  for (class in unique(classes)) {
    # classify of the all subtypes to the rest      
    is_class <- df[,classes == class]
    rest <- df[,classes != class]
    # calculate the mean of each samples in each of the genes
    mean_class <- data.frame(apply(is_class, 1, mean))
    mean_rest  <- data.frame(apply(rest, 1, mean))
    # rename the mean colunms
    colnames(mean_class) <- c('mean_class') 
    colnames(mean_rest) <- c('mean_rest')
    # range the expression of them
    mean_class <- arrange(mean_class,desc(mean_class))
    mean_rest <-arrange(mean_rest,desc(mean_rest))
    # change the probe names to genes
    #probe_to_gene(mean_class,'mean')
    #probe_to_gene(mean_rest,'mean')
    # save the files
    write.table(mean_class,file = sprintf('data/top_or_buttom_25/%s_class_GBM.txt', class))
    write.table(mean_rest,file = sprintf('data/top_or_buttom_25/%s_rest_GBM.txt',class))
  }
}

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
  {write(overlap, file =sprintf('data/top_or_buttom_25/0.45/%s_signatures_GBM.txt', percentage, subtypes))}}


percentage <- 0.45
# test
for(subtypes in unique(GBM_classes)){Overlap_genes(subtypes, percentage)}
