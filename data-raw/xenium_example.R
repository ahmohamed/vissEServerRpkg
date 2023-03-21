xenium_example = xenium(
  h5 = "data-raw/xenium/cell_feature_matrix.h5",
  tissue_positions = "data-raw/xenium/cells.csv",
  collections = c('CP:REACTOME', 'CP:BIOCARTA', 'CP:PID', 'CP:KEGG', 'CP:WIKIPATHWAYS', 'h', 'GO:BP', 'GO:CC', 'GO:MF'))

usethis::use_data(xenium_example, overwrite = TRUE)
jsonlite::write_json(jsonlite::fromJSON(xenium_example), '../examples/xenium.json')
