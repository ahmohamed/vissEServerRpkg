#----single-cell readers----
read_3files <- function(matrix.loc, gene.loc, barcode.loc) {
  tryCatch({
    gene.info <- read.delim(gene.loc, header = FALSE, colClasses = "character",
    stringsAsFactors = FALSE, quote = "", comment.char = "")
  }, error=function(x) {
      stop("Invalid features file. Features file should be a tab-delmited file with 3 columns (or its gz version)")
    }
  )
  if (ncol(gene.info) != 3) {
    stop("Invalid features file. Features file should be a tab-delmited file with 3 columns (or its gz version)")
  }
  colnames(gene.info) <- c("ID", "Symbol", "Type")

  tryCatch(
    { mat = as(Matrix::readMM(matrix.loc), "dgCMatrix") },
    error=function(x) {
      stop("Invalid matrix file.")
    }
  )
  barcodes = readLines(barcode.loc)
  bar_lengths = unique(sapply(barcodes, nchar))
  if (length(bar_lengths) > 1) {
    stop("Invalid barcode file. Barcodes should have the same length.")
  }

  if (bar_lengths > 50) {
    stop("Invalid barcode file. Barcodes should have be <50 bps.")
  }

  if ( !all(dim(mat) == c(nrow(gene.info), length(barcodes))) ) {
    stop("Matrix file does not match features and barcode files: Incompatible dimensions.")
  }

  list(
    mat = mat,
    cell.names = barcodes,
    gene.info = gene.info
  )
}

sc10x <- function (mat, features, barcodes, col.names = TRUE) {
  sample.names = list("Sample")
  # load.out <- DropletUtils:::.read_from_hdf5(h5, genome = genome, version = 'auto')
  load.out <- list(read_3files(mat, features, barcodes))

  nsets <- 1
  full_data <- vector("list", nsets)
  gene_info_list <- vector("list", nsets)
  cell_info_list <- vector("list", nsets)
  for (i in seq_len(nsets)) {
    current <- load.out[[i]]
    full_data[[i]] <- current$mat
    gene_info_list[[i]] <- current$gene.info
    cell.names <- current$cell.names
    cell_info_list[[i]] <- S4Vectors::DataFrame(Sample = rep(sample.names[i],
      length(cell.names)), Barcode = cell.names, row.names = NULL)
  }
  if (nsets > 1 && length(unique(gene_info_list)) != 1L) {
    stop("gene information differs between runs")
  }
  gene_info <- gene_info_list[[1]]
  rownames(gene_info) <- gene_info$ID
  full_data <- do.call(cbind, full_data)
  cell_info <- do.call(rbind, cell_info_list)
  if (col.names) {
    if (nsets == 1L) {
      cnames <- cell_info$Barcode
    }
    else {
      sid <- rep(seq_along(cell_info_list), vapply(cell_info_list,
        nrow, 1L))
      cnames <- paste0(sid, "_", cell_info$Barcode)
    }
    colnames(full_data) <- cnames
  }
  SingleCellExperiment::SingleCellExperiment(list(counts = full_data), rowData = gene_info,
    colData = cell_info, metadata = list(Samples = list("Sample")))
}

readSce <- function(mat, features, barcodes, annot = NULL){
  sce <- sc10x(mat, features, barcodes, col.names = TRUE)
  # change names to Symbol
  rownames(sce) <- SummarizedExperiment::rowData(sce)$Symbol
  sce
}

#' @import SpatialExperiment

#----spatial readers----
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
