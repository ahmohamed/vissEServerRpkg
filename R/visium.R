#----xenium/cosmx----
fun_if <- function(x, run, fun, ...) {
  if (run)
    fun(x, ...)
  else
    x
}

#' @export
#' @param expr file path to h5 (xenium), .*_exprMat_file.csv.gz file (CosMx) or
#'   matrix.mtx.gz (scRNA-seq).
#' @param tissue_positions file path to cells.csv (xenium) or
#'   .*_metadata_file.csv
cellular <- funwrapper(function(expr,
                                tissue_positions = NA_character_,
                                features = NA_character_,
                                barcodes = NA_character_,
                                tech = c('single-cell', 'visium', 'xenium', 'cosmx'),
                                method_filter = 'fixed',
                                method_normalise = 'none',
                                method_features = 'all',
                                sum = 10,
                                detected = 10,
                                neg = 5,
                                min_gene_count = 0,
                                dimred = c('PCA', 'NMF', 'UMAP', 'TSNE'),
                                dimred_fa = c('PCA', 'NMF'),
                                ncomponents = 5,
                                top_n_sets = 1000,
                                idtype = 'SYM',
                                org = 'hs',
                                collections = 'all') {

  set.seed(1234)
  tech = match.arg(tech)
  dimred_fa = match.arg(dimred_fa)

  # read files
  message('Reading files')
  if (tech == 'single-cell') {
    if (is.na(features) | is.na(barcodes))
      stop(sprintf('Files missing for %s data', tech))
    spe = readSce(expr, features, barcodes)
    neg_regex = '^MT|^mt'
  } else if (tech == 'visium') {
    if (is.na(tissue_positions))
      stop(sprintf('Files missing for %s data', tech))
    pos = readVisiumXY(tissue_positions)
    spe = readSpe(expr, pos)
    neg_regex = '^MT|^mt'
  } else if (tech == 'xenium') {
    if (is.na(tissue_positions))
      stop(sprintf('Files missing for %s data', tech))
    pos = readXeniumXY(tissue_positions)
    spe = readSpe(expr, pos)
    neg_regex = '^BLANK|^NegControl'
  } else if (tech == 'cosmx') {
    if (is.na(tissue_positions))
      stop(sprintf('Files missing for %s data', tech))
    spe = readCosMx(expr, tissue_positions)
    neg_regex = '^NegPrb'
  }

  #process data
  message('Preprocessing data')
  spe = spe |>
    addQC(neg_regex = neg_regex) |>
    filterCells(
      method = method_filter,
      sum = sum,
      detected = detected,
      neg = neg
    ) |>
    dropLowCount(lower = min_gene_count) |>
    dropNeg(neg_regex = neg_regex) |>
    doNormalise(method = method_normalise) |>
    doSelectFeatures(method = method_features) |>
    fun_if('PCA' %in% dimred, runPCA, ncomponents = 50) |>
    fun_if('NMF' %in% dimred, runNMF, ncomponents = 20) |>
    fun_if('UMAP' %in% dimred, runUMAP) |>
    fun_if('TSNE' %in% dimred, runTSNE)

  out = scVisseFA(
    sce = spe,
    dimred = dimred_fa,
    ncomponents = ncomponents,
    top_n_sets = top_n_sets,
    idtype = idtype,
    org = org,
    collections = collections
  )

  message('Serializing Results')
  out$api_version = api_version
  out$method = ifelse(tech == 'single-cell', 'SC', 'visium')
  out
})
