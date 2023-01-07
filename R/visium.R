#' @import SpatialExperiment

read_tissue_positions <- function(tissue_positions) {
  df <- read.csv(tissue_positions, header = FALSE, row.names = 1)
  if (ncol(df) != 5) {
    stop("Invalid tissue_positions file. ",
      "This should be a CSV file with 6 columns")
  }
  colnames(df) <- c("in_tissue", "array_row", "array_col",
    "pxl_col_in_fullres", "pxl_row_in_fullres")

  df$in_tissue <- as.logical(df$in_tissue)
  df
}

readSpe <- function(hdf5file, tissue_positions, name="Visium Sample") {
  tryCatch({
    sce <- DropletUtils::read10xCounts(samples = hdf5file, sample.names = name, col.names = TRUE, type="HDF5")
  }, error=function(x) {
    stop("Could not read H5 file.")
  })

  spd <- SpatialExperiment:::.read_xyz(tissue_positions)
  if (ncol(sce) != nrow(spd)) {
    warning("Data file does not match tissue position file: Incompatible dimensions.")
  }

  obs <- intersect(colnames(sce), rownames(spd))
  if (length(obs) < 100) {
    stop('There are less than 100 tissue positions with valid data.')
  }

  sce <- sce[, obs]
  spd <- spd[obs, ]
  spe = SpatialExperiment::SpatialExperiment(
    assays = SummarizedExperiment::assays(sce),
    rowData = S4Vectors::DataFrame(Symbol = SummarizedExperiment::rowData(sce)$Symbol),
    sample_id = name, spatialData = S4Vectors::DataFrame(spd),
    spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"))

  rownames(spe) <- SummarizedExperiment::rowData(spe)$Symbol
  spe
}

#' @export
visium <- funwrapper(function(h5, tissue_positions,
  filter_cell = "adaptive", sum=NULL, detected=NULL, mito=NULL,
  hvg=2000, min_gene_count = 0,
  dimred="PCA", ncomponents=5, top_n_sets=1000,
  idtype='SYM', org='hs', collections='all'
) {
  message("Reading files")
  spe_raw <- readSpe(h5, tissue_positions)
  xyranges <- list(xcoord=range(spatialCoords(spe_raw)[,1]), ycoord = range(spatialCoords(spe_raw)[,2]))

  message("Preprocessing data")
  # subset the object to keep only spots over tissue
  spe <- spe_raw[, spatialData(spe_raw)$in_tissue]
  spe <- scPreprocess(spe, filter_cell=filter_cell, sum=sum,
    detected=detected, mito=mito, hvg=hvg, min_gene_count=min_gene_count) %>%
    runPCA()

  message("Performing FA")
  msigdb = getCollections(idtype=idtype, org=org, collections=collections)
  out = visseFA(sce=spe, msigdb=msigdb, dimred=dimred, ncomponents=ncomponents, top_n_sets=top_n_sets, org=org)
  out = summarizeFA(out)
  # out = readRDS("server/R/robjects/examples/visiumFA.RDS")
  message("Serializing Results")
  out$api_version = api_version
  out$method = "visium"
  out
  # lapply(out, jsonlite::toJSON, digits=2)
  # jsonlite::toJSON(SpatialExperiment:::.read_xyz(tissue_positions), digits=2)
})
