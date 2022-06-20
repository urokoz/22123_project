#### Course 22123: Computational precision medicine
#### Project work
#### By Mathias Rahbek-Borre, Yu Liu and Joachim Breitenstein
#### 20/06/2022
#### Technical University of Denmark

# Calculate which genes are significantly different for each subtype
diff_class <- c()

for (class in unique(CIT_classes)) {
  # classify of the all subtypes to the rest      
  is_class <- CIT_full[,CIT_classes == class]
  rest <- CIT_full[, CIT_classes != class]
  
  mean_class <- data.frame(apply(is_class, 1, mean))
  mean_rest  <- data.frame(apply(rest, 1, mean))
  
  colnames(mean_class) <- c('mean_class') 
  colnames(mean_rest) <- c('mean_rest')
  # range the expression of them
  mean_class <- arrange(mean_class,desc(mean_class))
  mean_rest <-arrange(mean_rest,desc(mean_rest))
  
  mean_class_names <- rownames(mean_class)
  mean_rest_names <- rownames(mean_rest)
  
  match_idx <- match(mean_class_names, mean_rest_names)
  
  rank_vector <- c()
  dist_vector <- c()
  for (i in 1:length(match_idx)) {
    rank_vector[i] <- as.integer(i)
    dist_vector[i] <- match_idx[i] - as.integer(i)
  }
  
  mean_class$rank <- rank_vector
  mean_class$diff <- dist_vector
  
  diff_class[[class]] <- arrange(mean_class, desc(diff))
}

pred_perf <- function(dataset, signa, classes, class) {
  enrich <- gsva(dataset,
                 signa,
                 method="ssgsea",
                 ssgsea.norm = F,
                 verbose=F)
  return(ks.test(enrich[classes == class], enrich[classes != class])$statistic)
}

signatures <- c()
performances <- c()

for (class in c("normL", "lumA", "lumB", "lumC", "basL", "mApo")) {
  sig_probes <- rownames(diff_class[[class]])
  
  best_perf <- 0
  conseq_worse <- 0
  best_size <- NULL
  for (i in 1:length(sig_probes)) {
    signa <- list(sig_probes[1:i])
    
    perf <- pred_perf(CIT_full, signa, CIT_classes, class)
    print(sprintf("%s, %d:  %f", class, i, perf))
    performances <- rbind(performances, t(c(class, i, perf)))
    
    if (perf <= best_perf) {
      conseq_worse <- conseq_worse + 1
      if (conseq_worse == 25) {
        break
      }
    } else {
      conseq_worse <- 0
      best_perf <- perf
      best_size <- i
    }
  }
  signatures[[class]] <- list(sig_probes[1:best_size])
}

save(signatures, file = "data/rank_difference/CIT_probe_signatures.Rdata")

df <- data.frame(performances)
colnames(df) <- c("Subtype", "Sig_len", "Perf")
df$Perf <- as.numeric(as.character(df$Perf))
df$Sig_len <- as.numeric(as.character(df$Sig_len))

df %>% 
  ggplot(aes(x = Sig_len,
             y = Perf,
             color = Subtype)) +
  geom_line() + 
  labs(title = "Performance of variation in subtypes' signature length", x = "Number of genes in signatures", y = "Performance")



## testing and visualization of difference in gene expression for one subtype vs. rest for one gene
# class <- "basL"
# gene_number <- 1
#
# test_class <- data.frame(class = class, val = CIT_full[rownames(diff_class[[class]])[gene_number], CIT_classes == class])
# test_rest <- data.frame(class = "rest", val = CIT_full[rownames(diff_class[[class]])[gene_number], CIT_classes != class])
# 
# rbind(test_class, test_rest) %>% 
#   boxplot(class, val, "Expression for subtype", "Subtype", "Expression")


# function version of this code:
rank_diff_fnc <- function(dataset, classes) {
  diff_class <- c()
  for (class in unique(classes)) {
    # classify of the all subtypes to the rest      
    is_class <- dataset[,classes == class]
    rest <- dataset[, classes != class]
    
    mean_class <- data.frame(apply(is_class, 1, mean))
    mean_rest  <- data.frame(apply(rest, 1, mean))
    
    colnames(mean_class) <- c('mean_class') 
    colnames(mean_rest) <- c('mean_rest')
    # range the expression of them
    mean_class <- arrange(mean_class,desc(mean_class))
    mean_rest <-arrange(mean_rest,desc(mean_rest))
    
    mean_class_names <- rownames(mean_class)
    mean_rest_names <- rownames(mean_rest)
    
    match_idx <- match(mean_class_names, mean_rest_names)
    
    rank_vector <- c()
    for (i in 1:length(match_idx)) {
      rank_vector[i] <- as.integer(i)
      dist_vector[i] <- match_idx[i] - as.integer(i)
    }
    
    mean_class$rank <- rank_vector
    mean_class$diff <- dist_vector
    
    diff_class[[class]] <- rownames(arrange(mean_class, desc(diff)))
  }
  return(diff_class)
}


eval_signatures <- function(dataset, genes_of_interest, classes) {
  signatures <- c()
  performances <- c()
  
  for (class in unique(classes)) {
    sig_probes <- genes_of_interest[[class]]
    
    best_perf <- 0
    conseq_worse <- 0
    best_size <- NULL
    for (i in 1:length(sig_probes)) {
      signa <- list(sig_probes[1:i])
      
      perf <- pred_perf(dataset, signa, classes, class)
      print(sprintf("%s, %d:  %f", class, i, perf))
      performances <- rbind(performances, t(c(class, i, perf)))
      
      if (perf <= best_perf) {
        conseq_worse <- conseq_worse + 1
        if (conseq_worse == 25) {
          break
        }
      } else {
        conseq_worse <- 0
        best_perf <- perf
        best_size <- i
      }
    }
    signatures[[class]] <- list(sig_probes[1:best_size])
  }
  return(signatures)
}

diff_test <- rank_diff_fnc(CIT_full, CIT_classes)

test_signatures <- calc_signatures(CIT_full, diff_test, CIT_classes)


for (class in unique(CIT_classes)) {
  signa <- data.frame(signatures[[class]])
  colnames(signa) <- c("Probe.Set.ID")
  
  df_joined <- Bordet_annot %>%
    select(Probe.Set.ID, Gene.Symbol) %>% 
    mutate(Gene.Symbol = str_extract(Gene.Symbol, "^\\S+")) %>%
    inner_join(signa, by = "Probe.Set.ID") %>% 
    filter(Gene.Symbol != "---")
  
  file_name <- sprintf("data/signatures3_CIT_%s_genes.txt", class)
  write.table(df_joined$Gene.Symbol, file = file_name, quote = F, row.names = F, col.names = F)
}




