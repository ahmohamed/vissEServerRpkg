#' @export
gsea = funwrapper(function(genelist, idtype='SYM', org='hs', collections='all'){
  message(sprintf(
    "Starting GSEA with %d genes, ID type %s, organism %s and collections %s",
    length(genelist), idtype, org, paste(collections, collapse = ", ")
  ))
  genelist = setNames(as.numeric(sapply(genelist, "[[", 2)), sapply(genelist, "[[", 1))
  msigdb = getCollections(idtype=idtype, org=org, collections=collections)
  message(sprintf("Testing enrichement for %d genesets", length(msigdb)))

  pathways = setNames(lapply(msigdb, GSEABase::geneIds), lapply(msigdb, GSEABase::setName))
  res <- fgsea::fgsea(pathways, genelist)
  res_sig = res[res$padj < 0.05 & !is.na(res$padj) & !is.na(res$pathway), ]
  message(sprintf("Found %d signifacntly enriched genesets", nrow(res_sig)))

  siggs = msigdb[res_sig$pathway]
  gsfdr = setNames(res_sig$padj, res_sig$pathway)

  gset_attrs = dplyr::select(res_sig, ID=pathway, FDR=padj, Count=size, NES)
  gset_attrs$FDR = -log10(gset_attrs$FDR)

  out = visseWrapper(siggs, -log10(gsfdr), genelist, "logFC", gset_attrs = gset_attrs, org=org)
  out$geneSummary = geneSummary(msigdb, names(genelist))
  out$api_version = api_version
  out$method = "GSEA"
  out
})
