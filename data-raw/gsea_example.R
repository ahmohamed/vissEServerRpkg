genes_df = readRDS("data-raw/emt_allg.rds") |>
  dplyr::filter(gene_biotype %in% 'protein_coding') |>
  dplyr::select(gene_name, logFC)

readr::write_csv(genes_df, 'data-raw/gsea/genes_table.csv')
df = split(genes_df, 1:nrow(genes_df))

gsea_example = gsea(df, collections = c('CP:REACTOME', 'CP:BIOCARTA', 'CP:PID', 'CP:KEGG', 'CP:WIKIPATHWAYS', 'h', 'GO:BP', 'GO:CC', 'GO:MF'))
usethis::use_data(gsea_example, overwrite = TRUE)
jsonlite::write_json(jsonlite::fromJSON(gsea_example), '../examples/gsea.json')
