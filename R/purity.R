## use the estimate to calculate the enrichment score a stromal 
## content gene set and an immune infiltrate gene set in agiven sample

# download the package 

library('estimate')
library('tidyverse')
library("broom")
source("R/99_func_file.R")

CIT_gene_master <- probe_to_gene("CIT", "max")
CIT_gene <- CIT_gene_master[2:nrow(CIT_gene_master),]

names(CIT_gene)[names(CIT_gene) == "Gene.Symbol"] <- "NAME"

CIT_gene <- column_to_rownames(CIT_gene, var = "NAME")

write.table(CIT_gene, file = "data/_raw/CIT_samples.txt", sep = "\t", quote = F)

filterCommonGenes(input.f="data/_raw/CIT_samples.txt", output.f='data/_raw/CIT_genes.gct', id="GeneSymbol")

# calculate the score 
estimateScore(input.ds = 'data/_raw/CIT_genes.gct',
              output.ds = 'data/_raw/CIT_scores.gct',
              platform = c("affymetrix"))

# read and tidy the output file
CIT_gene <- rownames_to_column(CIT_gene, var = "NAME")

purity_expr <- CIT_gene %>% 
  pivot_longer(cols = !NAME, names_to = "Samples", values_to = "Expr") %>% 
  pivot_wider(names_from = NAME, values_from = Expr) 

purity_data <- read_tsv('data/_raw/CIT_scores.gct', skip = 2)

purity_data <- purity_data[4, -2] %>% 
  pivot_longer(!NAME, names_to = "Samples", values_to = "Purity") %>% 
  select(!NAME)

test <- data.frame(CIT_classes) %>% 
  rownames_to_column(var = "Samples") %>% 
  inner_join(purity_data, by = "Samples")


purity_expr_master <- purity_data %>% 
  inner_join(purity_expr, by = "Samples")


purity_expr_master %>% 
  ggplot(aes(x = Purity,
             y = AACS)) +
  geom_point()


lin_mdl <- purity_expr_master %>% 
  pivot_longer(!c(Samples, Purity), names_to = "Gene", values_to = "Expression") %>% 
  group_by(Gene) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(mu = map(data,
                  ~glm(Expression ~ Purity,
                       data = .x)))

purity_corrected_CIT <- lin_mdl %>%
  mutate(resi = map(mu, residuals),
         expected_expr = map(mu,
                             ~predict(., data.frame(Purity = 1)))) %>% 
  unnest(c(data, resi, expected_expr)) %>% 
  select(!mu) %>% 
  mutate(corr_expr = (resi + expected_expr)) %>% 
  select(Samples, Gene, corr_expr) %>% 
  pivot_wider(names_from = Samples, values_from = corr_expr)

save(purity_corrected_CIT, file = "data/CIT_purity_corrected.Rdata")
  

test %>% boxplot(CIT_classes, Purity, "Purity of all subtypes", "Subtype", "Purity")

