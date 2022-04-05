#' @importFrom fgsea fgsea

#' @export
gsea = function(genelist, idtype='SYM', org='hs', collections='all'){
  genelist = setNames(as.numeric(sapply(genelist, "[[", 2)), sapply(genelist, "[[", 1))
  msigdb = getCollections(idtype=idtype, org=org, collections=collections)
  print(length(msigdb))
  pathways = setNames(lapply(msigdb, GSEABase::geneIds), lapply(msigdb, GSEABase::setName))
  res <- fgsea(pathways, genelist)
  res_sig = res[res$padj < 0.05, ]
  siggs = msigdb[res_sig$pathway]
  gsfdr = setNames(res_sig$padj, res_sig$pathway)

  out = visseWrapper(siggs, -log10(gsfdr), genelist, "logFC")
  out$api_version = api_version
  out$method = "GSEA"
  jsonlite::toJSON(out, digits=2)
}
