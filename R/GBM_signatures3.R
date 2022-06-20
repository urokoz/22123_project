#### Course 22123: Computational precision medicine
#### Project work
#### By Mathias Rahbek-Borre, Yu Liu and Joachim Breitenstein
#### 20/06/2022
#### Technical University of Denmark


high_diff_ranks_GBM <- rank_diff_fnc(GBM_expr, GBM_classes)

high_diff_signatures_GBM <- calc_signatures(GBM_expr, high_diff_ranks_GBM, GBM_classes)


performances <- high_diff_signatures_GBM$performances
signatures <- high_diff_signatures_GBM$signatures

save(signatures, file = "data/rank_difference/GBM_high_diff_signatures.Rdata")

for (class in unique(GBM_classes)) {
  signa <- data.frame(signatures[[class]])
  colnames(signa) <- c("Gene.Symbol")
  
  file_name <- sprintf("data/rank_difference/signature_GBM_%s_genes.txt", class)
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
# 
# 
# # Calculate which genes are significantly different for each subtype
# diff_class <- c()
# for (class in unique(GBM_classes)) {
#   # classify of the all subtypes to the rest      
#   is_class <- GBM_expr[,GBM_classes == class]
#   rest <- GBM_expr[,GBM_classes != class]
#   
#   mean_class <- data.frame(apply(is_class, 1, mean))
#   mean_rest  <- data.frame(apply(rest, 1, mean))
#   
#   colnames(mean_class) <- c('mean_class') 
#   colnames(mean_rest) <- c('mean_rest')
#   # range the expression of them
#   mean_class <- arrange(mean_class,desc(mean_class))
#   mean_rest <-arrange(mean_rest,desc(mean_rest))
#   
#   mean_class_names <- rownames(mean_class)
#   mean_rest_names <- rownames(mean_rest)
#   
#   match_idx <- match(mean_class_names, mean_rest_names)
#   
#   rank_vector <- c()
#   for (i in 1:length(match_idx)) {
#     rank_vector[i] <- as.integer(i)
#     dist_vector[i] <- match_idx[i] - as.integer(i)
#   }
#   
#   mean_class$rank <- rank_vector
#   mean_class$diff <- dist_vector
#   
#   diff_class[[class]] <- arrange(mean_class, desc(diff))
# }
# 
# signatures <- c()
# performances <- c()
# 
# for (class in unique(GBM_classes)) {
#   sig_probes <- rownames(diff_class[[class]])
#   
#   best_perf <- 0
#   conseq_worse <- 0
#   best_size <- NULL
#   for (i in 1:length(sig_probes)) {
#     signa <- list(sig_probes[1:i])
#     
#     perf <- pred_perf(GBM_expr, signa, GBM_classes, class)
#     print(sprintf("%s, %d:  %f", class, i, perf))
#     performances <- rbind(performances, t(c(class, i, perf)))
#     
#     if (perf <= best_perf) {
#       conseq_worse <- conseq_worse + 1
#       if (conseq_worse == 25) {
#         break
#       }
#     } else {
#       conseq_worse <- 0
#       best_perf <- perf
#       best_size <- i
#     }
#   }
#   signatures[[class]] <- list(sig_probes[1:best_size])
# }
# 
# save(signatures, file = "data/rank_difference/GBM_gene_signatures.Rdata")
# 
# 
# rank_diff_fnc(GBM_expr, GBM_classes)
# 
# calc_signatures(GBM_expr, diff_test, GBM_classes)
# 
# 
