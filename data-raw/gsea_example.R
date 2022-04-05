library(tidyverse)
genes_df = readRDS("data-raw/emt_allg.rds") %>%
  filter(gene_biotype %in% 'protein_coding') %>%
  dplyr::select(gene_name, logFC)


df = split(genes_df, 1:nrow(genes_df))

gsea_example = gsea(df)
usethis::use_data(gsea_example, overwrite = TRUE)
