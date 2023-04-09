api_version = 0.1

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

getCollections <- function(idtype='SYM', org='hs', collections='all', minSize=0, maxSize=100000) {
  if (!idtype %in% c('SYM', 'EZ')) {
    idtype = 'SYM'
  }
  msigdb = getdata(paste0(org, '_', idtype))
  if (!'all' %in% collections) {
    msigdb = subsetCollection(msigdb, collections)
  }
  if (minSize > 0 || maxSize < 100000) {
    keep = sapply(msigdb, \(x) length(GSEABase::geneIds(x)) >= minSize && length(GSEABase::geneIds(x)) <= maxSize)
    msigdb = msigdb[keep]
  }
  msigdb
}

handle_ids <- function(ids, msigdb, org, idtype) {
  if (idtype %in% c('SYM', 'EZ')) {
    return (setNames(ids, ids))
  }

  sep_groups = data.frame(ids=ids) |>
    dplyr::mutate(original_ids=ids) |>
    tidyr::separate_rows(ids, sep = ';')

  orgdblist = list(hs=org.Hs.eg.db::org.Hs.eg.db, mm=org.Mm.eg.db::org.Mm.eg.db)
  symbol_ids = to_symbol(sep_groups$ids, orgdblist[[org]], from = idtype)
  universe = unique(unlist(GSEABase::geneIds(msigdb)))
  sep_groups$symbols = symbol_ids[sep_groups$ids]
  mapped_symbols = sep_groups |>
    dplyr::mutate(matched=!is.na(symbols) & symbols %in% universe) |>
    dplyr::group_by(original_ids) |> dplyr::arrange(!matched) |>
    dplyr::slice_head(n=1) |>
    dplyr::mutate(symbols=dplyr::coalesce(symbols, ids))

  setNames(mapped_symbols$symbols, mapped_symbols$original_ids)
}

to_symbol <- function(ids, orgdb, from) {
  if (!is(orgdb, 'OrgDb'))
    stop('Provide valid organism database')

  if (!any(ids %in% AnnotationDbi::keys(orgdb, from)))
    stop(sprintf('None of the "%s" IDs are valid', from))

  #convert numeric IDs to character
  if (!is.character(ids))
    ids = as.character(ids)

  #convert types
  ids = AnnotationDbi::mapIds(orgdb, ids, 'SYMBOL', from, multiVals = 'first')

  return(ids)
}

getPPI <- function(org="hs") {
  imex_0821 = getdata('imex_0821')
  taxid = c(hs="9606",  mm="10090")[org]
  ppi = imex_0821[imex_0821$Taxid %in% taxid, ]
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
    'edges' =  igraph::as_data_frame(gs_ovnet, 'edges'),
    'nodes' = dplyr::left_join(grps_df, vertices_df, by=c(Geneset = 'label')),
    # 'groups' = plyr::ldply(grps, function(x) data.frame('Geneset' = x), .id = 'Cluster'),
    # 'geneMembership' = head(geneIds(siggs[V(gs_ovnet)$label])),
    'words' = words
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
