
#' @export
ora = funwrapper(function(
  genelist,
  org,
  idtype,
  collections='all',
  minSize=0,
  maxSize=100000,
  pvalue=0.05,
  overlap_measure = c("ari", "jaccard", "ovlapcoef"),
  thresh=0.25) {
  message(sprintf(
    "Starting ORA with %d genes, ID type %s, organism %s and collections %s",
    length(genelist), idtype, org, paste(collections, collapse = ", ")
  ))
  msigdb = getCollections(org=org, collections=collections)
  message(sprintf("Testing enrichement for %d genesets", length(msigdb)))

  genelist = handle_ids(ids=unlist(genelist), gsc=msigdb, org=org, idtype=idtype)
  genelist = genelist[!duplicated(genelist)]

  gene_summary = geneSummary(msigdb, genelist)

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

  #prepare for cluster profiler
  msigdb_cp = lapply(msigdb, function(x)
    data.frame('term' = GSEABase::setName(x), 'gene' = GSEABase::geneIds(x)))
  msigdb_cp = do.call(rbind, msigdb_cp)
  res = clusterProfiler::enricher(genelist, TERM2GENE = msigdb_cp, minGSSize=minSize, maxGSSize=maxSize)
  if (is.null(res)) {
    stop("Pathway enrichment failed. Check the provided gene list.")
  }

  res = as.data.frame(res)
  res_sig = res[res$p.adjust < pvalue & !is.na(res$p.adjust), ]
  if (nrow(res_sig) < 10) {
    stop("Enrichment results has less that 10 significant genesets. Consider adding more geneset collections.")
  }

  message(sprintf("Found %d significantly enriched genesets", nrow(res_sig)))

  gset_attrs = dplyr::select(res_sig, ID, FDR=p.adjust, Count)
  gset_attrs$FDR = -log10(gset_attrs$FDR)

  siggs = msigdb[res_sig$ID]
  gsfdr = setNames(res$p.adjust, res_sig$ID)
  out = visseWrapper(siggs, -log10(gsfdr), gset_attrs = gset_attrs, org=org, overlap_measure=overlap_measure, thresh=thresh)
  out$gene_summary = gene_summary
  out$geneset_summary = genesetSummary(msigdb, out)
  out$api_version = api_version
  out$method = "ORA"
  out
})
