#' @import SpatialExperiment

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

#----visium----
#' @export
visium <- funwrapper(function(h5,
                              tissue_positions,
                              method_filter = 'adaptive',
                              method_normalise = 'scran',
                              method_features = 'HVG',
                              sum = 1000,
                              detected = 300,
                              mito = 10,
                              n_features = 2000,
                              min_gene_count = 0,
                              dimred = 'PCA',
                              ncomponents = 5,
                              top_n_sets = 1000,
                              idtype = 'SYM',
                              org = 'hs',
                              collections = 'all') {
  set.seed(1234)
  message('Reading files')
  pos = readVisiumXY(tissue_positions)
  spe = readSpe(h5, pos)
  # subset the object to keep only spots over tissue
  spe = spe[, spatialData(spe)$in_tissue]

  message('Preprocessing data')
  spe = spe |>
    addQC() |>
    filterCells(
      method = method_filter,
      sum = sum,
      detected = detected,
      neg = mito
    ) |>
    dropLowCount(min_gene_count) |>
    doNormalise(method = method_normalise) |>
    doSelectFeatures(method = method_features, n = n_features) |>
    runPCA(ncomponents = 50) |>
    runNMF(ncomponents = 20) |>
    runUMAP() |>
    runTSNE()

  out = scVisseFA(
    sce = spe,
    dimred = dimred,
    ncomponents = ncomponents,
    top_n_sets = top_n_sets,
    idtype = idtype,
    org = org,
    collections = collections
  )

  message('Serializing Results')
  out$api_version = api_version
  out$method = 'visium'
  out
})

#----xenium----
#' @export
xenium <- funwrapper(function(h5,
                              tissue_positions,
                              method_filter = 'fixed',
                              method_normalise = 'none',
                              method_features = 'all',
                              sum = 10,
                              detected = 10,
                              neg = 5,
                              neg_regex = '^BLANK|^NegControl',
                              min_gene_count = 0,
                              dimred = 'PCA',
                              ncomponents = 5,
                              top_n_sets = 1000,
                              idtype = 'SYM',
                              org = 'mm',
                              collections = 'all') {
  set.seed(1234)
  message('Reading files')
  pos = readXeniumXY(tissue_positions)
  spe = readSpe(h5, pos)

  message('Preprocessing data')
  spe = spe |>
    addQC(neg_regex = neg_regex) |>
    filterCells(
      method = method_filter,
      sum = sum,
      detected = detected,
      neg = neg
    ) |>
    dropLowCount(min_gene_count)

  #drop control probes
  spe = spe[!grepl(neg_regex, rownames(spe)), ]
  spe = spe |>
    doNormalise(method = method_normalise) |>
    doSelectFeatures(method = method_features) |>
    runPCA(ncomponents = 50) |>
    runNMF(ncomponents = 20) |>
    runUMAP() |>
    runTSNE()

  out = scVisseFA(
    sce = spe,
    dimred = dimred,
    ncomponents = ncomponents,
    top_n_sets = top_n_sets,
    idtype = idtype,
    org = org,
    collections = collections
  )

  message('Serializing Results')
  out$api_version = api_version
  out$method = 'xenium'
  out
})
