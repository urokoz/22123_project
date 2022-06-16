#### in this section I will calculate and rank the enrichment of different genes in
#### sample level and see if there are some genes are constantly on the top 25% of the
#### specific subtype and in the button of the rest.

## load the expression data
## the CIT_full contains all the expression data for all probes 
## the CIT_annot contains annotions of the 375 probes used for classification  
## the CIT class contains the subtypes of all the samples

load('data/_raw/CIT_data.Rdata')
load('data/_raw/Bordet.rdata')
source("R/99_func_file.R")
CIT_full <- probe_to_gene('CIT', 'mean')

# load some necessary packages
library('tidyverse')


# Calculate which genes are significantly different for each subtype

for (class in unique(CIT_classes)) {
# classify of the all subtypes to the rest      
    is_class <- CIT_full[,CIT_classes == class]
    rest <- CIT_full[, CIT_classes != class]
# calculate the mean of each sanples in each of the genes
    mean_class <- data.frame(apply(is_class, 1, mean))
    mean_rest  <- data.frame(apply(rest, 1, mean))
# rename the mean colunms
    colnames(mean_class) <- c('mean_class') 
    colnames(mean_rest) <- c('mean_rest')
# range the expression of them
    mean_class <- arrange(mean_class,desc(mean_class))
    mean_rest <-arrange(mean_rest,desc(mean_rest))
#save the files
    write.table(mean_class,file = sprintf('data/top_or_buttom_25/%s_class.txt', class))
    write.table(mean_rest,file = sprintf('data/top_or_buttom_25/%s_rest.txt',class))}



