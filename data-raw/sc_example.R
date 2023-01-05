sc_example = sc(
  mat = "data-raw/sc/matrix.mtx.gz",
  features = "data-raw/sc/features.tsv.gz",
  barcodes = "data-raw/sc/barcodes.tsv",
  collections = c('CP:REACTOME', 'CP:BIOCARTA', 'CP:PID', 'CP:KEGG', 'CP:WIKIPATHWAYS', 'h', 'GO:BP', 'GO:CC', 'GO:MF'))

usethis::use_data(sc_example, overwrite = TRUE)
jsonlite::write_json(jsonlite::fromJSON(sc_example), '../examples/sc.json')

