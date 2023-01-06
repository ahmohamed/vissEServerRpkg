scPreprocess <- function(sce,
  filter_cell = "adaptive", sum=NULL, detected=NULL, mito=NULL,
  hvg=2000, min_gene_count = 0
){
  sce %>%
    doQC() %>%
    filterCells(method = filter_cell, sum=sum, detected=detected, mito=mito) %>%
    dropLowCount(min_gene_count) %>%
    doNormalise() %>%
    doSelectFeatures(n=hvg)
}

doQC <- function(sce){
  # get mito genes
  is_mito = grepl("(^MT-)|(^mt-)", SummarizedExperiment::rowData(sce)$Symbol)
  if (!is.null(SummarizedExperiment::rowData(sce)$Chr)) {
    is_mito = is_mito | SummarizedExperiment::rowData(sce)$Chr %in% 'MT'
  }

  mito_genes <- SummarizedExperiment::rowData(sce)$Symbol[is_mito]
  scater::addPerCellQC(sce, subsets = list(Mito = mito_genes))
}

cellThresholds <- function(sce, method = "adaptive", sum=NULL, detected=NULL, mito=NULL){
  method = match.arg(method)
  if(! all(c("sum", "detected", "subsets_Mito_percent") %in% names(SummarizedExperiment::colData(sce)))){
    stop("sce does not have QC stats")
  }

  if(method == "adaptive"){
    scater::isOutlier(sce$sum, log = TRUE, type = "lower") |
      scater::isOutlier(sce$detected, log = TRUE, type = "lower") |
      scater::isOutlier(sce$subsets_Mito_percent, log = TRUE, type = "higher")
  } else {
    sce$sum < threshold$sum |
      sce$detected < threshold$detected |
      sce$subsets_Mito_percent > threshold$mito_percent
  }
}

filterCells <- function(sce,  method = "adaptive", sum=NULL, detected=NULL, mito=NULL){
  sce[, !cellThresholds(sce,  method, sum, detected, mito)]
}

doNormalise <- function(sce, ...){
  set.seed(1234)
  clust_sce <- scran::quickCluster(sce, min.size = min(100, ncol(sce)))
  sce <- scran::computeSumFactors(sce, cluster = clust_sce, min.mean=0.1)
  scater::logNormCounts(sce)
}

doSelectFeatures <- function(sce, n = 2000, ...){
  dec <-  scran::modelGeneVar(sce) # HVG 2000
  SingleCellExperiment::rowSubset(sce, "HVG") <-  scran::getTopHVGs(dec, n = n)

  sce
}

dropLowCount <- function(sce, lower = 0){
  sce[scuttle::perFeatureQCMetrics(sce)$detected > lower, ]
}

## Dimred
runPCA <- function(sce, ncomponents = 50){
  # Use HVGs if present
  genes <- SummarizedExperiment::rowData(sce)$Symbol
  if(!is.null(SummarizedExperiment::rowData(sce)$HVG)){
    genes <- genes[SummarizedExperiment::rowData(sce)$HVG]
  }

  scater::runPCA(sce, ncomponents = 50, subset_row = genes)
}

runUMAP <- function(sce){
  scater::runUMAP(sce, dimred = "PCA")
}



## FA
doFA <- function(sce, gsc, dimred="PCA", ncomponents=15) {
  pca_weights <- attr(SingleCellExperiment::reducedDim(sce, type = dimred), "rotation")[, 1:ncomponents]
  gsc = gsc[sapply(gsc, function(x) sum(GSEABase::geneIds(x) %in% rownames(pca_weights)) >= 5)]
  franks = singscore::rankGenes(pca_weights)
  pal_fsea = singscore::multiScore(franks, gsc)$Scores
  pal_fsea
}

visseFA <- function(sce, msigdb, dimred, ncomponents=15, top_n_sets=1000, cutoff_scores=0.3, org='hs') {
  if (top_n_sets < 100) {
    stop("top_n_sets should be at least 100 geneset")
  }

  gene_summary = geneSummary(msigdb, rownames(sce))
  if (gene_summary$value[["Used in Analysis"]] == 0) {
    stop("No genes were mapped. Provided genelist is invalid.")
  }

  message(sprintf(
    "%d out of %d genes provided will be used for the analysis with %d genes not mapped",
    gene_summary$value[["Used in Analysis"]], length(rownames(sce)), gene_summary$value[["Not Mapped"]]
  ))
  pal_fsea <- doFA(sce, msigdb, dimred=dimred, ncomponents=ncomponents)

  gene_counts = lapply(GSEABase::geneIds(msigdb), intersect, rownames(sce)) |> lapply(length)
  gset_attrs = data.frame(ID=names(gene_counts), Count=unlist(gene_counts))

  lapply(1:ncomponents, function(fct) {
    fct_weights = attr(SingleCellExperiment::reducedDim(sce, type = dimred), "rotation")[, fct]
    fct_proj = SingleCellExperiment::reducedDim(sce, dimred)[, fct]
    fct_sc = pal_fsea[, fct]
    top_n_sets = head(sort(abs(fct_sc), decreasing = TRUE), top_n_sets)
    siggs = msigdb[names(top_n_sets)]
    out = visseWrapper(siggs, gsStats=fct_sc, gStats=fct_weights, gStat_name="Weight", gset_attrs = gset_attrs, org=org)
    out$gene_summary = gene_summary
    out$geneset_summary = genesetSummary(msigdb, out)
    out
  })
}


summarizeFA <- function(fa_results) {
  out = list(results=fa_results)
  out$summary = lapply(fa_results, function(f) {
    f$words %>%
      dplyr::filter(as.numeric(Cluster) < 5) %>%
      dplyr::group_by(Cluster) %>%
      dplyr::slice_max(freq, n=20, with_ties = FALSE)
  })
  out
}
