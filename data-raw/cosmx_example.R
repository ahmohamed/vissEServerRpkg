cosmx_example = cellular(.silent=FALSE, 
  expr = "data-raw/cosmx/Lung5_Rep3_exprMat_file.csv.gz",
  tissue_positions = "data-raw/cosmx/Lung5_Rep3_metadata_file.csv",
  tech = 'cosmx',
  dimred = c('PCA', 'UMAP'),
  dimred_fa = 'PCA',
  collections = c('CP:REACTOME', 'CP:BIOCARTA', 'CP:PID', 'CP:KEGG', 'CP:WIKIPATHWAYS', 'h', 'GO:BP', 'GO:CC', 'GO:MF'),
  org = 'hs')

usethis::use_data(cosmx_example, overwrite = TRUE)
jsonlite::write_json(jsonlite::fromJSON(cosmx_example), '../examples/cosmx.json')
