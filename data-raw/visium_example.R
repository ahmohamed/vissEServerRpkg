library(tidyverse)
visium_example = visium(
  h5 = "data-raw/visium/raw_feature_bc_matrix.h5",
  tissue_positions = "data-raw/visium/tissue_positions_list.csv",
  collections = c('CP:REACTOME', 'CP:BIOCARTA', 'CP:PID', 'CP:KEGG', 'CP:WIKIPATHWAYS', 'h', 'GO:BP', 'GO:CC', 'GO:MF'))

usethis::use_data(visium_example, overwrite = TRUE)
