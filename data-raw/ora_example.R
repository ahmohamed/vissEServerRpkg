genes = readRDS("data-raw/emt_allg.rds") |>
  dplyr::filter(gene_biotype %in% 'protein_coding', FDR < 0.01)

readr::write_csv(genes, 'data-raw/ora/genes_table.csv')
ora_example = ora(genes$gene_name, collections = c('CP:REACTOME', 'CP:BIOCARTA', 'CP:PID', 'CP:KEGG', 'CP:WIKIPATHWAYS', 'h', 'GO:BP', 'GO:CC', 'GO:MF'))
usethis::use_data(ora_example, overwrite = TRUE)
jsonlite::write_json(jsonlite::fromJSON(ora_example), '../examples/ora.json')
