api_version = 0.1

getSpecies <- function() {
  fpath = system.file("extdata/species_info.rds", package = "vissEServer")
  species = readRDS(fpath)
  species$Species
}

getIdTypes <- function(org = getSpecies()) {
  org = match.arg(org)

  #load supported ID data
  fpath = system.file("extdata/species_geneid_present.rds", package = "vissEServer")
  idtypes = readRDS(fpath)

  # subset organism
  if (org %in% c('hsapiens', 'mmusculus')) {
    idtypes = colnames(idtypes)
  } else {
    idtypes = idtypes[org, ]
    idtypes = names(idtypes)[idtypes]
  }
  idtypes
}

getdata <- function(x) {
  e <- new.env()
  name <- utils::data(list=c(x), envir = e)[[1]]
  e[[name]]
}

subsetCollection <- function(gsc, collection = c()) {
  stopifnot(length(gsc) > 0)
  #filter collection & sub-collection
  ctype = lapply(gsc, GSEABase::collectionType)
  gsc = gsc[sapply(ctype, GSEABase::bcCategory) %in% collection |
      sapply(ctype, GSEABase::bcSubCategory) %in% collection]

  return(gsc)
}

getCollections <- function(org = getSpecies(), collections = "all", minSize = 0, maxSize = 100000) {
  org = match.arg(org)

  # load gene-set collection
  fpath = system.file(sprintf("extdata/species_gsc/%s_gsc.rds", org), package = "vissEServer")
  gsc = readRDS(fpath)

  # subset collections
  if (!"all" %in% collections) {
    gsc = subsetCollection(gsc, collections)
  }

  # subset by size
  if (minSize > 0 || maxSize < 100000) {
    keep = sapply(gsc, \(x) length(GSEABase::geneIds(x)) >= minSize && length(GSEABase::geneIds(x)) <= maxSize)
    gsc = gsc[keep]
  }
  gsc
}

handle_ids <- function(ids, gsc, org = getSpecies(), idtype = getIdTypes(org)) {
  org = match.arg(org)
  idtype = match.arg(idtype)

  # split gene groups (proteomics)
  sep_groups = data.frame(ids = ids) |>
    dplyr::mutate(original_ids = ids) |>
    tidyr::separate_rows(ids, sep = ";")

  # if idtype does not match database key, convert
  if (org %in% c("hsapiens", "mmusculus") & idtype != "external_gene_name") {
    sep_groups$mapped_ids = convert_ids(ids, org, idtype, "external_gene_name")
  } else if (!org %in% c("hsapiens", "mmusculus") & idtype != "ensembl_gene_id") {
    sep_groups$mapped_ids = convert_ids(ids, org, idtype, "ensembl_gene_id")
  } else {
    # if idtype matches key type, do nothing
    sep_groups$mapped_ids = sep_groups$ids
  }
  
  #define universe using gsc
  universe = unique(unlist(GSEABase::geneIds(gsc)))
  sep_groups = sep_groups |>
    dplyr::mutate(matched = !is.na(mapped_ids) & mapped_ids %in% universe) |>
    dplyr::group_by(original_ids) |>
    dplyr::arrange(!matched) |>
    dplyr::slice_head(n = 1) |>
    dplyr::mutate(mapped_ids = dplyr::coalesce(mapped_ids, ids))

  setNames(sep_groups$mapped_ids, sep_groups$original_ids)
}

convert_ids <- function(ids, org = getSpecies(), from = getIdTypes(org), to = getIdTypes(org)) {
  org = match.arg(org)
  from = match.arg(from)
  to = match.arg(to)

  # load gene ID map for organism
  fpath = system.file(sprintf("extdata/species_gsc/%s_idmap.rds", org), package = "vissEServer")
  idmap = readRDS(fpath)

  # check keys
  if (!any(ids %in% idmap[, from])) {
    stop(sprintf('None of the "%s" IDs are valid', from))
  }

  # convert numeric IDs to character
  if (!is.character(ids)) {
    ids = as.character(ids)
  }

  # select first mapping and crete mapping vec
  idmap = idmap[, c(from, to)]
  idmap = idmap[!duplicated(idmap[, from]), ]
  idmap = setNames(idmap[, to], idmap[, from])
  #convert IDs
  ids = idmap[ids]

  return(ids)
}

qmax <- function(x, prob = 0.99) {
  x_max = stats::quantile(abs(x), probs = prob)
  x = pmin(x, x_max)
  x = pmax(x, -x_max)
  return(x)
}

nice <- function(x, prob = 0.99, digits = 3) {
  signif(qmax(x, prob = prob), digits = digits)
}

emptyPPI <- function() {
  data.frame(
    EntrezA = character(),
    EntrezB = character(),
    Taxid = character(),
    InteractorA = character(),
    InteractorB = character(),
    SymbolA = character(),
    SymbolB = character(),
    InteractionType = character(),
    DetectionMethod = character(),
    Confidence = numeric(),
    Inferred = logical()
  )
}

getPPI <- function(org = getSpecies()) {
  org = match.arg(org)

  # taxid map
  taxid = c(hsapiens = "9606", mmusculus = "10090")
  if (!org %in% names(taxid)) {
    return(emptyPPI())
  } else {
    taxid = taxid[org]
  }

  #read PPI network
  fpath = system.file("extdata/imex.rds", package = "vissEServer")
  imex = readRDS(fpath)
  
  #subset to species
  ppi = imex[imex$Taxid %in% taxid, ]
  ppi
}

geneSummary <- function(msigdb, genes) {
  universe = unique(unlist(GSEABase::geneIds(msigdb)))
  both = sum(universe %in% genes)
  not_used = sum(!genes %in% universe)
  stats = list("Universe Size"=length(universe), "Used in Analysis"=both, "Not Mapped"=not_used)
  tibble::tibble(stat=names(stats), value=stats)
}

genesetSummary <- function(msigdb, out) {
  nodes = dplyr::ungroup(out$nodes)
  geneset_stats = nodes |> dplyr::summarise(
    "Significant Genesets"=dplyr::n(), "Categories"=dplyr::n_distinct(Category),
    "Subcategories"=dplyr::n_distinct(SubCategory), "Average Geneset Size"=mean(Size, na.rm=T),
    "Geneset Clusters"=dplyr::n_distinct(Cluster)
  ) |> unlist() |> floor()
  avg_cluster_size = nodes |> dplyr::group_by(Cluster) |> dplyr::count() |> dplyr::pull(n) |> mean(na.rm=T) |> floor()

  stats = as.list(c(
    "Tested Genesets"=length(msigdb),
    geneset_stats,
    "Average Cluster Size"=avg_cluster_size
  ))

  list(
    stats=tibble::tibble(stat=names(stats), value=stats),
    category_tally = nodes |> dplyr::group_by(Category) |> dplyr::summarise(count=dplyr::n()),
    subcategory_tally = nodes |> dplyr::group_by(Category, SubCategory) |> dplyr::summarise(count=dplyr::n())
  )
}

visseWrapper <- function(siggs, gsStats, gStats = NULL, gStat_name="Gene-level statistic", gset_attrs=NULL, org="hs",
  thresh=0.25, overlap_measure = c("ari", "jaccard", "ovlapcoef")) {
  message(sprintf("Starting vissE analysis for %d genesets", length(siggs)))
  if (length(siggs) < 10) {
    stop("Enrichment results has less that 10 significant genesets. Consider adding more geneset collections.")
  }

  #compute geneset overlaps between significant genesets
  gs_ovlap = vissE::computeMsigOverlap(siggs, thresh = thresh, measure=overlap_measure)
  message(sprintf("Detected %d genesets overlaps with threshold %f", nrow(gs_ovlap), thresh))

  #create a network from overlaps
  gs_ovnet = vissE::computeMsigNetwork(gs_ovlap, siggs)
  message(sprintf("Geneset network computed with %d nodes and %d edges", igraph::vcount(gs_ovnet), igraph::ecount(gs_ovnet)))

  #identify clusters
  grps = igraph::cluster_walktrap(gs_ovnet)
  #extract clustering results
  grps = igraph::groups(grps)
  message(sprintf("Found %d geneset clusters", length(grps)))

  #sort by stat
  sortstat = plyr::ldply(grps, function(x) {
    fdr.stat = median(abs(gsStats[x]))
    size.stat = length(x)
    return(c(fdr.stat, size.stat))
  })[, -1]
  sortstat = apply(apply(sortstat, 2, rank), 1, prod)
  grps = grps[order(sortstat, decreasing = TRUE)]
  names(grps) = 1:length(grps)
  grps_df = purrr::map_dfr(grps, function(x) data.frame('Geneset' = x), .id = 'Cluster')

  # use numeric ids instead of names to reduce output size
  igraph::V(gs_ovnet)$label = igraph::V(gs_ovnet)$name
  igraph::V(gs_ovnet)$name = 1:igraph::vcount(gs_ovnet)
  igraph::V(gs_ovnet)$degree = igraph::degree(gs_ovnet)
  vertices_df = igraph::as_data_frame(gs_ovnet, 'vertices') |>
    dplyr::mutate(SubCategory = dplyr::coalesce(SubCategory, Category))
  if (!is.null(gset_attrs)) {
    vertices_df = dplyr::left_join(vertices_df, gset_attrs, by=c("label"="ID"))
  }

  message("Computing wordclouds")
  p2 = suppressWarnings(vissE::plotMsigWordcloud(siggs, grps, type = 'Name', org = org))
  words = dplyr::rename(p2$data[, 1:3], Cluster = NodeGroup)
  words$Cluster = gsub(' \\(.*\\)', '', as.character(words$Cluster))
  out = list(
    'edges' =  igraph::as_data_frame(gs_ovnet, 'edges') |> dplyr::mutate(weight=signif(weight, 3)),
    'nodes' = dplyr::left_join(grps_df, vertices_df, by=c(Geneset = 'label')) |> dplyr::mutate_if(is.double, signif, digits=3),
    # 'groups' = plyr::ldply(grps, function(x) data.frame('Geneset' = x), .id = 'Cluster'),
    # 'geneMembership' = head(geneIds(siggs[V(gs_ovnet)$label])),
    'words' = dplyr::mutate(words, freq=signif(freq, 3))
  )
  if (!is.null(gStats)) {
    ppi = getPPI(org=org)

    #compute gene-level stats
    message(sprintf("Computing PPI network for %d genes", length(gStats)))
    p3 = vissE::plotGeneStats(gStats, siggs, grps, statName = gStat_name, topN = 5)
    out$genestats = p3$data |> dplyr::arrange(Group, rank) |> dplyr::select(-rank)

    comp_ppi = vissE:::computeMsigGroupPPI(ppi, siggs, grps, gStats, org=org)
    if (igraph::ecount(comp_ppi) > 0) {
      ppi_grps = comp_ppi |>
        tidygraph::activate(nodes) |>
          tidygraph::filter(!is.na(Degree)) |> tidygraph::select(group=Group, val=Degree, name=label) |>
        tidygraph::activate(edges) |> tidygraph::select(from, to, inferred=Inferred) |>
        tidygraph::to_split(group, split_by = 'nodes') |> lapply(as.list) |> setNames(NULL)

      out$ppi_grps = ppi_grps
    }
  }

  out
}
