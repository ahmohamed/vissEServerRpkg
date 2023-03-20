#' @export
gsea = funwrapper(function(genelist, idtype='SYM', org='hs', collections='all', scoretype="std", minsize=3){
  message(sprintf(
    "Starting GSEA with %d genes, ID type %s, organism %s and collections %s",
    length(genelist), idtype, org, paste(collections, collapse = ", ")
  ))
  msigdb = getCollections(idtype=idtype, org=org, collections=collections)
  message(sprintf("Testing enrichement for %d genesets", length(msigdb)))

  genelist = setNames(as.numeric(sapply(genelist, "[[", 2)), sapply(genelist, "[[", 1))
  converted_ids = handle_ids(ids=names(genelist), msigdb=msigdb, org=org, idtype=idtype)
  names(genelist) = converted_ids[names(genelist)]
  genelist = genelist[!duplicated(names(genelist))]
  gene_summary = geneSummary(msigdb, names(genelist))

  message(sprintf(
    "%d out of %d genes provided will be used for the analysis with %d genes not mapped",
    gene_summary$value[["Used in Analysis"]], length(genelist), gene_summary$value[["Not Mapped"]]
  ))

  if (gene_summary$value[["Not Mapped"]] > length(genelist) * 0.9) {
    stop("More than 90% of genes were not mapped. ",
      "Please check the provided data matches the selected ID type and organism.")
  }

  if (gene_summary$value[["Not Mapped"]] > length(genelist) * 0.5) {
    warning("More than 50% of genes were not mapped. ",
      "You may need to check the provided data matches the selected ID type and organism.")
  }

  pathways = setNames(lapply(msigdb, GSEABase::geneIds), lapply(msigdb, GSEABase::setName))
  res <- fgsea::fgsea(pathways, genelist, nproc=1, scoreType=scoretype, minSize=minsize)
  res_sig = res[res$padj < 0.05 & !is.na(res$padj) & !is.na(res$pathway), ]
  message(sprintf("Found %d significantly enriched genesets", nrow(res_sig)))

  siggs = msigdb[res_sig$pathway]
  gsfdr = setNames(res_sig$padj, res_sig$pathway)

  gset_attrs = dplyr::select(res_sig, ID=pathway, FDR=padj, Count=size, NES)
  gset_attrs$FDR = -log10(gset_attrs$FDR)

  out = visseWrapper(siggs, -log10(gsfdr), genelist, "logFC", gset_attrs = gset_attrs, org=org)
  out$gene_summary = gene_summary
  out$geneset_summary = genesetSummary(msigdb, out)
  out$api_version = api_version
  out$method = "GSEA"
  out
})
