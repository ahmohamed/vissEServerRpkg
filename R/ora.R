
#' @export
ora = funwrapper(function(genelist, idtype='SYM', org='hs', collections='all') {
  message(sprintf(
    "Starting ORA with %d genes, ID type %s, organism %s and collections %s",
    length(genelist), idtype, org, paste(collections, collapse = ", ")
  ))
  msigdb = getCollections(idtype=idtype, org=org, collections=collections)
  message(sprintf("Testing enrichement for %d genesets", length(msigdb)))

  #prepare for cluster profiler
  msigdb_cp = lapply(msigdb, function(x)
    data.frame('term' = GSEABase::setName(x), 'gene' = GSEABase::geneIds(x)))
  msigdb_cp = do.call(rbind, msigdb_cp)
  res = clusterProfiler::enricher(genelist, TERM2GENE = msigdb_cp)
  res = as.data.frame(res)
  res_sig = res[res$p.adjust < 0.05 & !is.na(res$p.adjust), ]
  message(sprintf("Found %d significantly enriched genesets", nrow(res_sig)))

  gset_attrs = dplyr::select(res_sig, ID, FDR=p.adjust, Count)
  gset_attrs$FDR = -log10(gset_attrs$FDR)

  siggs = msigdb[res_sig$ID]
  gsfdr = setNames(res$p.adjust, res_sig$ID)
  out = visseWrapper(siggs, -log10(gsfdr), gset_attrs = gset_attrs, org=org)
  out$geneSummary = geneSummary(msigdb, genelist)
  out$api_version = api_version
  out$method = "ORA"
  out
})
