library(tidyverse)
genes = readRDS("data-raw/emt_allg.rds") %>%
  filter(gene_biotype %in% 'protein_coding') %>%
  slice_max(n=1000, order_by = logFC) %>% pull(gene_name)

ora_example = ora(genes, collections = c('CP:REACTOME', 'CP:BIOCARTA', 'CP:PID', 'CP:KEGG', 'CP:WIKIPATHWAYS', 'h', 'GO:BP', 'GO:CC', 'GO:MF'))
usethis::use_data(ora_example, overwrite = TRUE)
