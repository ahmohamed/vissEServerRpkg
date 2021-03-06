
#' @export
ora = function(genelist, idtype='SYM', org='hs', collections='all') {
  msigdb = getCollections(idtype=idtype, org=org, collections=collections)

  #prepare for cluster profiler
  msigdb_cp = lapply(msigdb, function(x)
    data.frame('term' = GSEABase::setName(x), 'gene' = GSEABase::geneIds(x)))
  msigdb_cp = do.call(rbind, msigdb_cp)
  res = clusterProfiler::enricher(genelist, TERM2GENE = msigdb_cp)
  res = as.data.frame(res)
  res_sig = res[res$p.adjust < 0.05, ]
  gset_attrs = dplyr::select(res_sig, ID, FDR=p.adjust, Count)
  gset_attrs$FDR = -log10(gset_attrs$FDR)

  siggs = msigdb[res_sig$ID]
  gsfdr = setNames(res$p.adjust, res_sig$ID)
  out = visseWrapper(siggs, -log10(gsfdr), gset_attrs = gset_attrs)
  out$api_version = api_version
  out$method = "ORA"
  jsonlite::toJSON(out, digits=2)
}
