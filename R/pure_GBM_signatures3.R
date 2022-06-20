#### Course 22123: Computational precision medicine
#### Project work
#### By Mathias Rahbek-Borre, Yu Liu and Joachim Breitenstein
#### 20/06/2022
#### Technical University of Denmark


high_diff_ranks_GBM <- rank_diff_fnc(GBM_expr, GBM_classes)

high_diff_signatures_GBM <- calc_signatures(GBM_expr, high_diff_ranks_GBM, GBM_classes)


performances <- high_diff_signatures_GBM$performances
signatures <- high_diff_signatures_GBM$signatures

save(signatures, file = "data/rank_difference/pure/GBM_high_diff_signatures.Rdata")

for (class in unique(GBM_classes)) {
  signa <- data.frame(signatures[[class]])
  colnames(signa) <- c("Gene.Symbol")
  
  file_name <- sprintf("data/rank_difference/pure/signature_GBM_%s_genes.txt", class)
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
