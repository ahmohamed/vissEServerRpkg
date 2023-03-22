xenium_example = cellular(
  expr = "data-raw/xenium/cell_feature_matrix.h5",
  tissue_positions = "data-raw/xenium/cells.csv",
  tech = 'xenium',
  dimred = c('PCA', 'UMAP'),
  dimred_fa = 'PCA',
  collections = c('CP:REACTOME', 'CP:BIOCARTA', 'CP:PID', 'CP:KEGG', 'CP:WIKIPATHWAYS', 'h', 'GO:BP', 'GO:CC', 'GO:MF'),
  org = 'mm')

usethis::use_data(xenium_example, overwrite = TRUE)
jsonlite::write_json(jsonlite::fromJSON(xenium_example), '../examples/xenium.json')
