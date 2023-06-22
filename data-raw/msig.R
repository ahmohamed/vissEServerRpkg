library(usethis)
library(msigdb)

#msigdb
hsapiens = appendKEGG(getMsigdb('hs'))
mmusculus = appendKEGG(getMsigdb("mm"))

saveRDS(hsapiens, './inst/extdata/species_gsc/hsapiens_gsc.rds')
saveRDS(mmusculus, "./inst/extdata/species_gsc/mmusculus_gsc.rds")
