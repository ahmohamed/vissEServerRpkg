#' @import SpatialExperiment

#----reading functions----
readVisiumXY <- function(tissue_positions) {
  df <- read.csv(tissue_positions, header = FALSE, row.names = 1)
  if (ncol(df) != 5) {
    stop("Invalid tissue_positions file. ",
         "This should be a CSV file with 6 columns")
  }
  colnames(df) <- c("in_tissue", "array_row", "array_col",
                    "y", "x")

  df$in_tissue <- as.logical(df$in_tissue)

  return(df)
}

#read positions for xenium
readXeniumXY <- function(tissue_positions) {
  spd = read.csv(tissue_positions, header = TRUE, row.names = 1)
  colnames(spd)[1:2] = c('x', 'y')

  return(spd)
}

#read positions for cosmx
readCosMx <- function(exprmat, metadata) {
  annots = read.csv(metadata, header = TRUE)
  colnames(annots)[7:8] = c('x', 'y')
  emat = read.csv(exprmat, header = TRUE)
  emat = emat[emat$cell_ID > 0, ]

  #add unique cell IDs
  rownames(emat) = paste(emat$fov, emat$cell_ID, sep = '_')
  rownames(annots) = paste(annots$fov, annots$cell_ID, sep = '_')
  emat = t(emat[, -(1:2)])

  if(!all(colnames(emat) %in% rownames(annots))) {
    stop('Annotations do not match the data.')
  }

  if (ncol(emat) != nrow(annots)) {
    warning('Data file does not match tissue position file: Incompatible dimensions.')
  }

  obs = intersect(colnames(emat), rownames(annots))
  if (length(obs) < 100) {
    stop('There are less than 100 tissue positions with valid data.')
  }

  emat = emat[, obs]
  annots = annots[obs, ]
  spe = SpatialExperiment::SpatialExperiment(
    assays = list('counts' = emat),
    colData = as(annots, 'DataFrame'),
    rowData = S4Vectors::DataFrame(Symbol = rownames(emat)),
    sample_id = 'Spatial Sample',
    spatialCoordsNames = c('x', 'y')
  )

  rownames(spe) = SummarizedExperiment::rowData(spe)$Symbol
  spe
}

readSpe <- function(hdf5file, positions) {
  tryCatch({
    sce = DropletUtils::read10xCounts(samples = hdf5file, sample.names = 'Spatial Sample', col.names = TRUE, type='HDF5')
  }, error=function(x) {
    stop('Could not read H5 file.')
  })

  if (ncol(sce) != nrow(positions)) {
    warning('Data file does not match tissue position file: Incompatible dimensions.')
  }

  obs = intersect(colnames(sce), rownames(positions))
  if (length(obs) < 100) {
    stop('There are less than 100 tissue positions with valid data.')
  }

  sce = sce[, obs]
  positions = positions[obs, ]
  spe = SpatialExperiment::SpatialExperiment(
    assays = SummarizedExperiment::assays(sce),
    rowData = S4Vectors::DataFrame(Symbol = SummarizedExperiment::rowData(sce)$Symbol),
    sample_id = 'Spatial Sample',
    spatialCoords = as.matrix(positions[, c('x', 'y')])
  )

  rownames(spe) = SummarizedExperiment::rowData(spe)$Symbol
  spe
}

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
