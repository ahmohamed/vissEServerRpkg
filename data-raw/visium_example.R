visium_example = cellular(.silent=FALSE,
  expr = "data-raw/visium/raw_feature_bc_matrix.h5",
  tissue_positions = "data-raw/visium/tissue_positions_list.csv",
  tech = 'visium',
  dimred = c('PCA', 'UMAP', 'TSNE'),
  dimred_fa = 'PCA',
  collections = c('CP:REACTOME', 'CP:BIOCARTA', 'CP:PID', 'CP:KEGG', 'CP:WIKIPATHWAYS', 'h', 'GO:BP', 'GO:CC', 'GO:MF'))

usethis::use_data(visium_example, overwrite = TRUE)
jsonlite::write_json(jsonlite::fromJSON(visium_example), '../examples/visium.json')
