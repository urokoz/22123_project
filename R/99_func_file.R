library(tidyverse)

##function for conversion of probe ID to gene ID 
#input: 1. argument: string of dataframe; either "bordet" or "CIT", 2. argument: string with method (max, mean or median) 
#output: table containing genes and desired calculation of expression values

load_function <- function(df) {
  df <- 
  return(df)
}


probe_to_gene <- function(data, method) {
  
  if (data == "bordet") {
    Bordet_array_df <- data.frame("Probe.Set.ID" = row.names(Bordet_array), Bordet_array)
    Bordet_df_joined <- Bordet_annot %>%
      mutate(Gene.Symbol = str_extract(Gene.Symbol, "^\\S+")) %>% 
      select(Probe.Set.ID, Gene.Symbol) %>% 
      inner_join(Bordet_array_df, by = "Probe.Set.ID")
    
    if (method == "max") {
      final_table <- Bordet_df_joined %>% 
        select(!Probe.Set.ID) %>% 
        pivot_longer(cols = !Gene.Symbol, names_to = "Sample", values_to = "Expression") %>% 
        group_by(Gene.Symbol, Sample) %>% 
        summarise(Max = max(Expression)) %>% 
        pivot_wider(names_from = Sample, values_from = Max)
    }
    
    if (method == "mean") {
      final_table <- Bordet_df_joined %>%  
        select(!Probe.Set.ID) %>% 
        pivot_longer(cols = !Gene.Symbol, names_to = "Sample", values_to = "Expression") %>% 
        group_by(Gene.Symbol, Sample) %>% 
        summarise(Mean = mean(Expression)) %>% 
        pivot_wider(names_from = Sample, values_from = Mean)
    }
    
    if (method == "median") {
      final_table <- Bordet_df_joined %>%  
        select(!Probe.Set.ID) %>% 
        pivot_longer(cols = !Gene.Symbol, names_to = "Sample", values_to = "Expression") %>% 
        group_by(Gene.Symbol, Sample) %>% 
        summarise(Median = median(Expression)) %>% 
        pivot_wider(names_from = Sample, values_from = Median)
    }
  }  

  if (data == "CIT") {
    CIT_array_df <- data.frame("Probe.Set.ID" = row.names(CIT_full), CIT_full)
    CIT_df_joined <- Bordet_annot %>%
      mutate(Gene.Symbol = str_extract(Gene.Symbol, "^\\S+")) %>%
      select(Probe.Set.ID, Gene.Symbol) %>% 
      inner_join(CIT_array_df, by = "Probe.Set.ID")
    
    if (method == "max") {
      final_table <- CIT_df_joined %>% 
        select(!Probe.Set.ID) %>% 
        pivot_longer(cols = !Gene.Symbol, names_to = "Sample", values_to = "Expression") %>% 
        group_by(Gene.Symbol, Sample) %>% 
        summarise(Max = max(Expression)) %>% 
        pivot_wider(names_from = Sample, values_from = Max)
    }
    
    if (method == "mean") {
      final_table <- CIT_df_joined %>% 
        select(!Probe.Set.ID) %>% 
        pivot_longer(cols = !Gene.Symbol, names_to = "Sample", values_to = "Expression") %>% 
        group_by(Gene.Symbol, Sample) %>% 
        summarise(Mean = mean(Expression)) %>% 
        pivot_wider(names_from = Sample, values_from = Mean)
    }
    
    if (method == "median") {
      final_table <- CIT_df_joined %>% 
        select(!Probe.Set.ID) %>% 
        pivot_longer(cols = !Gene.Symbol, names_to = "Sample", values_to = "Expression") %>% 
        group_by(Gene.Symbol, Sample) %>% 
        summarise(Median = median(Expression)) %>% 
        pivot_wider(names_from = Sample, values_from = Median)
      
    }
    
  }
  return(final_table[2:nrow(final_table),])
}


signature_procsess <- function(df, specification, method) {
  if (specification == "CIT") {
    df_joined <- Bordet_annot %>%
      select(Probe.Set.ID, Gene.Symbol) %>% 
      inner_join(df, by = "Probe.Set.ID")
  
    if (method == "max") {
      final_table <- df_joined %>% 
        group_by(Gene.Symbol) %>% 
        summarise(Max_expr = max(last_col())) 
    }
    
    if (method == "mean") {
      final_table <- df_joined %>% 
        group_by(Gene.Symbol) %>% 
        summarise(Mean_expr = mean(last_col())) 
    }
    
    if (method == "Median") {
      final_table <- df_joined %>%
        group_by(Gene.Symbol) %>% 
        summarise(Median_expr = median(last_col())) 
    }
  }
  
  return(final_table)
}


boxplot <- function(df, variable_dis, variable_con, title, xlab, ylab) {
  
  plt <- ggplot(df, mapping = aes(y = {{variable_con}},
                                  x = {{variable_dis}}, 
                                  fill = {{variable_dis}})) + 
    geom_violin() + 
    geom_boxplot(width = 0.05,
                 color = "black",
                 fill = "white",
                 outlier.shape = NA) +
    labs(title = title, 
         x = xlab,
         y = ylab) +
    theme_classic() +
    theme(legend.position="none")
  
  return(plt)
}

probeID_to_geneID <- function(df) {
  
  colnames(df) <- c("Probe.Set.ID", "Expression") 
  geneID <- df %>% left_join(CIT_annot, by = "Probe.Set.ID")
  
  return(geneID)
}

# df <- load()
# read.delim("data/top_or_bottom_25/normL_rest.txt")
# probeID_to_geneID()

significant_genes <- function(expr_data, classes) {
  
  signif_genes <- c()
  unique_classes <- unique(classes)
  for (class in unique(classes)) {
    
    results <- c()
    for (probe in rownames(expr_data)) {
      is_class <- expr_data[probe, classes == class]
      rest <- expr_data[probe, classes != class]
      res <- wilcox.test(x = as.numeric(is_class),
                         y = as.numeric(rest))
      results[[probe]] <- res$p.value
      
    }
    results <- p.adjust(results)
    results <- results[results < 0.05]
    
    signif_genes[[class]] <- results
    
  }
  return(signif_genes)
}

# save(significant_genes(CIT_full, CIT_classes), file = "data/CIT_signif_subtype_genes.Rdata")
# save(significant_genes(GBM_expr, GBM_clinical$GeneExp_Subtype), file = "data/GBM_signif_subtype_genes.Rdata")


FC_calc <- function(expr_data, signif_genes, classes) {
  
  FC_lists <- c()
  for (class in unique(classes)) {
    results <- signif_genes[[class]]
    
    FC_list <- c()
    
    for (probe in names(results)) {
      subtype_median <- median(as.numeric(expr_data[probe, classes == class]))
      rest_median <- median(as.numeric(expr_data[probe, classes != class]))
      
      FC_list <- rbind(FC_list, c(probe, log2(subtype_median/rest_median)))
      
    }
    colnames(FC_list) <- c("probe", "FC")
    FC_list <- data.frame(FC_list)
    signa <- arrange(data.frame(FC_list), desc(FC))
    
    file_name <- sprintf("data/ranked_CIT_%s_probes.txt", class)
    
    write.table(signa$probe, file = file_name, quote = F, row.names = F, col.names = F)
    
    FC_lists[[class]] <- signa$probes
    
  }
  return(FC_lists)
}


pred_perf <- function(data, signa, classes, class) {
  enrich <- gsva(data,
                 signa,
                 method="ssgsea",
                 ssgsea.norm = F,
                 verbose=F)
  return(ks.test(enrich[classes == class], enrich[classes != class])$statistic)
}


calc_signatures <- function(data, interest_genes_list, classes) {
  
  signatures <- c()
  performances <- c()
  
  for (class in unique(classes)) {
    interest_genes <- interest_genes_list[[class]]
    
    best_perf <- 0
    conseq_worse <- 0
    best_size <- NULL
    for (i in 1:length(interest_genes)) {
      signa <- list(interest_genes[1:i])
      
      perf <- pred_perf(data, signa, classes, class)
      print(sprintf("%s, %d:  %f", class, i, perf))
      performances <- rbind(performances, t(c(class, i, perf)))
      
      if (perf <= best_perf) {
        conseq_worse <- conseq_worse + 1
        if (conseq_worse == 50) {
          break
        }
      } else {
        conseq_worse <- 0
        best_perf <- perf
        best_size <- i
      }
    }
    signatures[[class]] <- list(interest_genes[1:best_size])
  }
  return(signatures)
}


# For each sample, the signature with the highest enrichment corresponds to the subtype you assign to the given sample
name_max <- function(column) {
  subtype <- names(column)[which.max(column)]
  return(subtype)
}

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
    write.table(mean_class,file = sprintf('data/top_or_buttom_25/%s_class.txt', class))
    write.table(mean_rest,file = sprintf('data/top_or_buttom_25/%s_rest.txt',class))
  }
}
