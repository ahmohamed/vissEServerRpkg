sc_example = cellular(.silent=FALSE, 
  expr = "data-raw/sc/matrix.mtx.gz",
  features = "data-raw/sc/features.tsv.gz",
  barcodes = "data-raw/sc/barcodes.tsv",
  tech = 'single-cell',
  method_filter ='adaptive',
  method_normalise = 'scran',
  method_features = 'HVG',
  dimred = c('PCA', 'UMAP'),
  dimred_fa = 'PCA',
  collections = c('CP:REACTOME', 'CP:BIOCARTA', 'CP:PID', 'CP:KEGG', 'CP:WIKIPATHWAYS', 'h', 'GO:BP', 'GO:CC', 'GO:MF'))

usethis::use_data(sc_example, overwrite = TRUE)
jsonlite::write_json(jsonlite::fromJSON(sc_example), '../examples/sc.json')

