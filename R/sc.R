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
  sample.names <- list("Sample")
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

  #add gene annotations
  AnnotationHub::setAnnotationHubOption("ASK", FALSE)
  ah <- AnnotationHub::AnnotationHub()

  # TODO: mm annotations
  ensdb  <-  ah[['AH89426']] # Ensembl 103 EnsDb for Homo sapiens
  SummarizedExperiment::rowData(sce)$Chr  <-  AnnotationDbi::mapIds(ensdb, rownames(sce), keytype = 'GENEID', column = 'SEQNAME')
  SummarizedExperiment::rowData(sce)$Biotype  <-  AnnotationDbi::mapIds(ensdb, rownames(sce), keytype = 'GENEID', column = 'GENEBIOTYPE')

  # change names to Symbol
  rownames(sce) <- SummarizedExperiment::rowData(sce)$Symbol
  sce
}

#' @export
sc <- funwrapper(function(mat, features, barcodes,
  filter_cell = "adaptive", sum=NULL, detected=NULL, mito=NULL,
  hvg=2000, min_gene_count = 0,
  dimred="PCA", ncomponents=5, top_n_sets=1000,
  idtype='SYM', org='hs', collections='all'
) {
  message("Reading files")
  sce <- readSce(mat, features, barcodes)

  message("Preprocessing data")
  # subset the object to keep only spots over tissue
  sce <- scPreprocess(sce, filter_cell=filter_cell, sum=sum,
    detected=detected, mito=mito, hvg=hvg, min_gene_count=min_gene_count) %>%
    runPCA()

  message("Performing FA")
  msigdb = getCollections(idtype=idtype, org=org, collections=collections)
  out = visseFA(sce=sce, msigdb=msigdb, dimred=dimred, ncomponents=ncomponents, top_n_sets=top_n_sets, org=org)
  out = summarizeFA(out)
  # out = readRDS("server/R/robjects/examples/visiumFA.RDS")
  message("Serializing Results")
  out$api_version = api_version
  out$method = "SC"
  out
  # lapply(out, jsonlite::toJSON, digits=2)
  # jsonlite::toJSON(SpatialExperiment:::.read_xyz(tissue_positions), digits=2)
})
