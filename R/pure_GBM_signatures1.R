#### Course 22123: Computational precision medicine
#### Project work
#### By Mathias Rahbek-Borre, Yu Liu and Joachim Breitenstein
#### 20/06/2022
#### Technical University of Denmark

signif_genes_GBM <- significant_genes(GBM_expr, GBM_classes)
save(signif_genes_GBM, file = "data/mann_whitney_u_test/pure/GBM_signif_subtype_genes.Rdata")

FC_GBM <- FC_calc(GBM_expr, signif_genes_GBM, GBM_classes)
save(FC_GBM, file = "data/mann_whitney_u_test/pure/GBM_FC.Rdata")

signatures_GBM <- calc_signatures(GBM_expr, FC_GBM, GBM_classes)

performances <- signatures_GBM$performances
signatures <- signatures_GBM$signatures

save(signatures, file = "data/mann_whitney_u_test/pure/GBM_signatures.Rdata")


for (class in unique(GBM_classes)) {
  signa <- data.frame(signatures[[class]])
  if (length(signa) < 1) {
    next
  }
  
  colnames(signa) <- c("Gene.Symbol")

  file_name <- sprintf("data/mann_whitney_u_test/pure/signature_GBM_%s_genes.txt", class)
  write.table(signa, file = file_name, quote = F, row.names = F, col.names = F)
}


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

