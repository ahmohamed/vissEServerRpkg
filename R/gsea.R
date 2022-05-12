#' @export
gsea = function(genelist, idtype='SYM', org='hs', collections='all'){
  genelist = setNames(as.numeric(sapply(genelist, "[[", 2)), sapply(genelist, "[[", 1))
  msigdb = getCollections(idtype=idtype, org=org, collections=collections)
  print(length(msigdb))
  pathways = setNames(lapply(msigdb, GSEABase::geneIds), lapply(msigdb, GSEABase::setName))
  res <- fgsea::fgsea(pathways, genelist)
  res_sig = res[res$padj < 0.05, ]
  siggs = msigdb[res_sig$pathway]
  gsfdr = setNames(res_sig$padj, res_sig$pathway)

  gset_attrs = dplyr::select(res_sig, ID=pathway, FDR=padj, Count=size, NES)
  gset_attrs$FDR = -log10(gset_attrs$FDR)

  out = visseWrapper(siggs, -log10(gsfdr), genelist, "logFC", gset_attrs = gset_attrs)
  out$api_version = api_version
  out$method = "GSEA"
  jsonlite::toJSON(out, digits=2)
}
