library(tidyverse)
library(parallel)
library(hash)

##function for conversion of probe ID to gene ID 
#input: 1. argument: string of dataframe; either "bordet" or "CIT", 2. argument: string with method (max, mean or median) 
#output: table containing genes and desired calculation of expression values

probe_to_gene <- function(data, method) {
  
  if (data == "bordet") {
    Bordet_array_df <- data.frame("Probe.Set.ID" = row.names(Bordet_array), Bordet_array)
    Bordet_df_joined <- Bordet_annot %>% left_join(Bordet_array_df, by = "Probe.Set.ID")
    
    if (method == "max") {
      final_table <- Brodet_df_joined %>% select(Gene.Symbol, HER2.13:ncol(Brodet_df_joined)) %>% 
        pivot_longer(cols = !Gene.Symbol, names_to = "Sample", values_to = "Expression") %>% 
        group_by(Gene.Symbol, Sample) %>% 
        summarise(Max = max(Expression)) %>% 
        pivot_wider(names_from = Sample, values_from = Max)
    }
    
    if (method == "mean") {
      final_table <- Brodet_df_joined %>% select(Gene.Symbol, HER2.13:ncol(Brodet_df_joined)) %>% 
        pivot_longer(cols = !Gene.Symbol, names_to = "Sample", values_to = "Expression") %>% 
        group_by(Gene.Symbol, Sample) %>% 
        summarise(Mean = mean(Expression)) %>% 
        pivot_wider(names_from = Sample, values_from = Mean)
    }
    
    if (method == "median") {
      final_table <- Brodet_df_joined %>% select(Gene.Symbol, HER2.13:ncol(Brodet_df_joined)) %>% 
        pivot_longer(cols = !Gene.Symbol, names_to = "Sample", values_to = "Expression") %>% 
        group_by(Gene.Symbol, Sample) %>% 
        summarise(Median = median(Expression)) %>% 
        pivot_wider(names_from = Sample, values_from = Median)
    }
  }  

  if (data == "CIT") {
    CIT_array_df <- data.frame("Probe.Set.ID" = row.names(CIT_full), CIT_full)
    CIT_df_joined <- Bordet_annot %>% left_join(CIT_array_df, by = "Probe.Set.ID")
    
    if (method == "max") {
      final_table <- CIT_df_joined %>% select(Gene.Symbol, CIT_DSOA_440:ncol(CIT_df_joined)) %>% 
        pivot_longer(cols = !Gene.Symbol, names_to = "Sample", values_to = "Expression") %>% 
        group_by(Gene.Symbol, Sample) %>% 
        summarise(Max = max(Expression)) %>% 
        pivot_wider(names_from = Sample, values_from = Max)
    }
    
    if (method == "mean") {
      final_table <- CIT_df_joined %>% select(Gene.Symbol, CIT_DSOA_440:ncol(CIT_df_joined)) %>% 
        pivot_longer(cols = !Gene.Symbol, names_to = "Sample", values_to = "Expression") %>% 
        group_by(Gene.Symbol, Sample) %>% 
        summarise(Mean = mean(Expression)) %>% 
        pivot_wider(names_from = Sample, values_from = Mean)
    }
    
    if (method == "median") {
      final_table <- CIT_df_joined %>% select(Gene.Symbol, CIT_DSOA_440:ncol(CIT_df_joined)) %>% 
        pivot_longer(cols = !Gene.Symbol, names_to = "Sample", values_to = "Expression") %>% 
        group_by(Gene.Symbol, Sample) %>% 
        summarise(Median = median(Expression)) %>% 
        pivot_wider(names_from = Sample, values_from = Median)
      
    }
    
  }
  return(final_table)
}

#test
#Brodet_df_joined %>% select(Gene.Symbol, HER2.13:ncol(Brodet_df_joined)) %>% 
 # pivot_longer(cols = !Gene.Symbol, names_to = "Sample", values_to = "Expression") %>% 
  #group_by(Gene.Symbol, Sample) %>% 
  #summarise(Median = median(Expression)) %>% 
  #pivot_wider(names_from = Sample, values_from = Median)



signature_procsess <- function(df, specification, method) {
  if (specification == "CIT") {
    df_joined <- inner_join(df, CIT_annot, by = "Probe.Set.ID")
  
    if (method == "max") {
      final_table <- df_joined %>% select(Gene.Symbol, 
                                          last_col()) %>% 
        group_by(Gene.Symbol) %>% 
        summarise(Max = max(last_col())) 
    }
    
    if (method == "mean") {
      final_table <- df_joined %>% select(Gene.Symbol, 
                                          last_col()) %>% 
        group_by(Gene.Symbol) %>% 
        summarise(Max = mean(last_col())) 
    }
    
    if (method == "Median") {
      final_table <- df_joined %>% select(Gene.Symbol, 
                                          last_col()) %>% 
        group_by(Gene.Symbol) %>% 
        summarise(Max = median(last_col())) 
    }
  }
  
  return(final_table)
}



