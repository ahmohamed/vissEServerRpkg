library(usethis)
for (version in getMsigdbVersions()[1:4]) {
  for (org in c('hs', 'mm')) {
    for (idtype in c('SYM', 'EZ')) {
      name = paste0(org, "_", idtype, "_", version)
      msigdb = getMsigdb(org = org, id = idtype, version=version)
      msigdb = appendKEGG(msigdb)
      assign(name, msigdb)

      # Weird behavior of use_data
      # See here https://stackoverflow.com/a/49676445/3192855
      do.call("use_data", list(as.name(name), overwrite = TRUE))
    }
  }
}


