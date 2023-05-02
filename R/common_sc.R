scPreprocess <- function(sce,
                         filter_cell = 'adaptive',
                         sum = NULL,
                         detected = NULL,
                         neg = NULL,
                         hvg = 2000,
                         min_gene_count = 0) {
  sce |>
    addQC() |>
    filterCells(
      method = filter_cell,
      sum = sum,
      detected = detected,
      neg = neg
    ) |>
    dropLowCount(min_gene_count) |>
    doNormalise() |>
    doSelectFeatures(n = hvg)
}

#----QC----
addQC <- function(sce, neg_regex = '(^MT-)|(^mt-)') {
  # get negative control genes (e.g., neg probes and mito)
  is_neg = grepl(neg_regex, SummarizedExperiment::rowData(sce)$Symbol)

  #use Chr if present
  if (!is.null(SummarizedExperiment::rowData(sce)$Chr)) {
    is_neg = is_neg | SummarizedExperiment::rowData(sce)$Chr %in% c('MT', 'mt')
  }

  neg_genes = SummarizedExperiment::rowData(sce)$Symbol[is_neg]
  scater::addPerCellQC(sce, subsets = list(Neg = neg_genes))
}

cellThresholds <-
  function(sce,
           method = c('adaptive', 'fixed'),
           sum = 0,
           detected = 0,
           neg = 0) {
    method = match.arg(method)
    if (!all(
      c('sum', 'detected', 'subsets_Neg_percent') %in% names(SummarizedExperiment::colData(sce))
    )) {
      stop('sce does not have QC stats')
    }

    if (method == 'adaptive') {
      scater::isOutlier(sce$sum, log = TRUE, type = 'lower') |
        scater::isOutlier(sce$detected, log = TRUE, type = 'lower') |
        scater::isOutlier(sce$subsets_Neg_percent, log = TRUE, type = 'higher')
    } else if (method == 'fixed') {
      sce$sum < sum |
        sce$detected < detected |
        sce$subsets_Neg_percent > neg
    } else {
      stop('threshold method not supported')
    }
  }

filterCells <-
  function(sce,
           method = 'adaptive',
           sum = NULL,
           detected = NULL,
           neg = NULL) {
    sce[, !cellThresholds(sce,  method, sum, detected, neg)]
  }

dropLowCount <- function(sce, lower = 0){
  sce[scuttle::perFeatureQCMetrics(sce)$detected > lower, ]
}

dropNeg <- function(sce, neg_regex = '(^MT-)|(^mt-)') {
  # get negative control genes (e.g., neg probes and mito)
  is_neg = grepl(neg_regex, SummarizedExperiment::rowData(sce)$Symbol)

  #use Chr if present
  if (!is.null(SummarizedExperiment::rowData(sce)$Chr)) {
    is_neg = is_neg | SummarizedExperiment::rowData(sce)$Chr %in% c('MT', 'mt')
  }

  sce[!is_neg, ]
}

#----normalisation----
normalise_none <- function(sce) {
  SingleCellExperiment::logcounts(sce) = log2(SingleCellExperiment::counts(sce) + 1)

  return(sce)
}

normalise_scran <- function(sce) {
  cl = scran::quickCluster(sce, min.size = min(100, ncol(sce)))
  sce = scran::computeSumFactors(sce, clusters = cl)
  sce = scater::logNormCounts(sce)

  return(sce)
}

normalise_sct <- function(sce) {
  vst_res = sctransform::vst(SingleCellExperiment::counts(sce, vst.flavor = 'v2'),
                             verbosity = 0)
  #subset data
  sce = sce[rownames(vst_res$y), ]
  SingleCellExperiment::logcounts(sce) = vst_res$y

  return(sce)
}

doNormalise <- function(sce, method = c('scran', 'none', 'sctransform'), ...){
  switch (match.arg(method),
          scran = normalise_scran(sce, ...),
          none = normalise_none(sce, ...),
          sctransform = normalise_sct(sce, ...)
  )
}

#----feature selection----
feature_HVGs <- function(sce, n = 2000) {
  stopifnot(is(sce, 'SingleCellExperiment'))

  #scran standard
  dec = scran::modelGeneVar(sce)
  features = scran::getTopHVGs(dec, n = n)
  SingleCellExperiment::rowSubset(sce, 'UseGenes') = features

  return(sce)
}

feature_all <- function(sce) {
  SingleCellExperiment::rowSubset(sce, 'UseGenes') = rownames(sce)

  return(sce)
}

doSelectFeatures <- function(sce, method = c('HVG', 'all'), ...){
  switch (
    match.arg(method),
    HVG = feature_HVGs(sce, ...),
    all = feature_all(sce, ...)
  )
}

#----dimred----
runPCA <- function(sce, ncomponents = 50) {
  scater::runPCA(
    sce,
    ncomponents = ncomponents,
    subset_row = SummarizedExperiment::rowData(sce)$UseGenes,
    name = 'PCA',
    exprs_values = 'logcounts'
  )
}

runNMF <- function(sce, ncomponents = 5) {
  sce = scater::runNMF(
    sce,
    ncomponents = ncomponents,
    subset_row = SummarizedExperiment::rowData(sce)$UseGenes,
    name = 'NMF',
    exprs_values = 'logcounts'
  )

  #add names
  colnames(SingleCellExperiment::reducedDim(sce, 'NMF')) = paste0('NMF', seq(ncomponents))
  colnames(attr(SingleCellExperiment::reducedDim(sce, 'NMF'), 'basis')) = paste0('NMF', seq(ncomponents))
  rownames(attr(SingleCellExperiment::reducedDim(sce, 'NMF'), 'basis')) = rownames(sce)[SummarizedExperiment::rowData(sce)$UseGenes]

  return(sce)
}

runUMAP <- function(sce) {
  scater::runUMAP(sce, dimred = 'PCA', name = 'UMAP')
}

runTSNE <- function(sce) {
  scater::runTSNE(sce, dimred = 'PCA', name = 'TSNE')
}

exportReducedDim = function(sce, name, ncomponents) {
  df = SingleCellExperiment::reducedDim(sce, name)
  colnames(df) = paste0(name, 1:ncol(df))
  if (ncomponents > ncol(df)) {
    ncomponents = ncol(df)
  }
  df[, 1:ncomponents] |> apply(2, nice)
}

#----clustering----
cluster_graph <- function(sce, snn = 20, alg = igraph::cluster_walktrap, ...) {
  if ('Cluster' %in% colnames(SummarizedExperiment::colData(sce))) {
    stop('Annotation column cannot be named "Cluster"')
  }

  g = scran::buildSNNGraph(sce, k = snn, use.dimred = 'PCA')
  cl = alg(g, ...)$membership
  sce$Cluster = as.factor(cl)

  return(sce)
}

#----factor interpretation----
doFA <- function(sce, gsc, dimred = c('PCA', 'NMF'), ncomponents = 15) {
  dimred = match.arg(dimred)
  directional = dimred %in% c('PCA') #are gene weights directional?
  attr_name = switch(
    dimred,
    PCA = 'rotation',
    NMF = 'basis',
    stop(sprintf('Factor analysis is not supported with %s', dimred))
  )

  gene_weights = attr(SingleCellExperiment::reducedDim(sce, type = dimred),
                      attr_name)[, 1:ncomponents]
  gsc = gsc[sapply(gsc, function(x)
    sum(GSEABase::geneIds(x) %in% rownames(gene_weights)) >= 5)]
  franks = singscore::rankGenes(gene_weights)
  fsea = singscore::multiScore(franks, gsc, centerScore = directional)$Scores
  attr(fsea, 'weights') = gene_weights
  fsea
}

visseFA <-
  function(sce,
           msigdb,
           dimred,
           ncomponents = 15,
           top_n_sets = 1000,
           cutoff_scores = 0.3,
           org = 'hs',
           thresh=0.25,
           overlap_measure = c("ari", "jaccard", "ovlapcoef")) {

    if (top_n_sets < 100) {
      stop('top_n_sets should be at least 100 geneset')
    }

    gene_summary = geneSummary(msigdb, rownames(sce))

    fsea <- doFA(sce, msigdb, dimred = dimred, ncomponents = ncomponents)

    gene_counts = lapply(GSEABase::geneIds(msigdb), intersect, rownames(sce)) |> lapply(length)
    gset_attrs = data.frame(ID=names(gene_counts), Count=unlist(gene_counts))

    lapply(1:ncomponents, function(fct) {
      fct_weights = attr(fsea, 'weights')[, fct]
      fct_proj = SingleCellExperiment::reducedDim(sce, dimred)[, fct]
      fct_sc = fsea[, fct]
      gset_attrs$Score = fct_sc[gset_attrs$ID]
      top_n_sets = head(sort(abs(fct_sc), decreasing = TRUE), top_n_sets)
      siggs = msigdb[names(top_n_sets)]
      out = visseWrapper(siggs, gsStats=signif(fct_sc, 3), gStats=signif(fct_weights, 3), gStat_name='Weight', gset_attrs = gset_attrs, org=org, overlap_measure=overlap_measure, thresh=thresh)
      out$gene_summary = gene_summary
      out$geneset_summary = genesetSummary(msigdb, out)
      out
    })
  }

summarizeFA <- function(fa_results) {
  out = list(results=fa_results)
  out$summary = lapply(fa_results, function(f) {
    f$words |>
      dplyr::filter(as.numeric(Cluster) < 5) |>
      dplyr::group_by(Cluster) |>
      dplyr::slice_max(freq, n=20, with_ties = FALSE)
  })
  out$gptsummary = lapply(fa_results, function(f) f$gpt)
  out
}

scVisseFA = function(sce,
                     filter_cell,
                     sum,
                     detected,
                     mito,
                     hvg,
                     min_gene_count,
                     dimred,
                     ncomponents,
                     top_n_sets,
                     org,
                     msigdb,
                     thresh=0.25,
                     overlap_measure = c("ari", "jaccard", "ovlapcoef")) {
  message('Performing factor interpretation')
  out = visseFA(
    sce = sce,
    msigdb = msigdb,
    dimred = dimred,
    ncomponents = ncomponents,
    top_n_sets = top_n_sets,
    org = org,
    overlap_measure=overlap_measure,
    thresh=thresh
  )
  out = summarizeFA(out)

  dimlist = lapply(SingleCellExperiment::reducedDimNames(sce), function(x)
    exportReducedDim(sce, name = x, ncomponents = ncomponents))
  if (is(sce, 'SpatialExperiment')) {
    coords = as.data.frame(spatialCoords(sce))
    colnames(coords) = c('Tissue X', 'Tissue Y')
    dimlist = append(dimlist, coords)
  }

  qcmetrics = SummarizedExperiment::colData(sce)[, c('sum', 'detected', 'subsets_Neg_percent')]
  cellmetrics = as.data.frame(cbind(qcmetrics, dimlist)) |>
    dplyr::mutate(
      sum = qmax(sum),
      detected = qmax(detected),
      subsets_Neg_percent = nice(subsets_Neg_percent)
    ) |>
    dplyr::rename(
      `Library Size` = sum,
      `Gene Count` = detected,
      `% Control/Mito Genes` = subsets_Neg_percent
    )
  rownames(cellmetrics) = NULL
  out$cellmetrics = as.list(cellmetrics)
  out
}
