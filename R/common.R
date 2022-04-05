#' @import vissE
#' @import msigdb
#' @import igraph
#' @import tidygraph
#' @import jsonlite

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

getCollections <- function(idtype='SYM', org='hs', collections='all') {
  msigdb = getdata(paste0(org, '_', idtype))
  if (collections != 'all') {
    msigdb = subsetCollection(msigdb, collections)
  }
  msigdb
}

getPPI <- function(org="hs") {
  imex_0821 = getdata('imex_0821')
  ppi = imex_0821[imex_0821$Taxid %in% '9606', ]
}

visseWrapper <- function(siggs, gsStats, gStats = NULL, gStat_name="Gene-level statistic") {
  #compute geneset overlaps between significant genesets
  gs_ovlap = computeMsigOverlap(siggs, thresh = 0.25)
  #create a network from overlaps
  gs_ovnet = computeMsigNetwork(gs_ovlap, siggs)

  #identify clusters
  grps = cluster_walktrap(gs_ovnet)
  #extract clustering results
  grps = igraph::groups(grps)
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
  V(gs_ovnet)$label = V(gs_ovnet)$name
  V(gs_ovnet)$name = 1:vcount(gs_ovnet)

  p2 = suppressWarnings(plotMsigWordcloud(siggs, grps, type = 'Name'))
  words = dplyr::rename(p2$data[, 1:3], Cluster = NodeGroup)
  words$Cluster = gsub(' \\(.*\\)', '', as.character(words$Cluster))
  out = list(
    'edges' =  igraph::as_data_frame(gs_ovnet, 'edges'),
    'nodes' = dplyr::left_join(grps_df, igraph::as_data_frame(gs_ovnet, 'vertices'), by=c(Geneset = 'label')),
    # 'groups' = plyr::ldply(grps, function(x) data.frame('Geneset' = x), .id = 'Cluster'),
    # 'geneMembership' = head(geneIds(siggs[V(gs_ovnet)$label])),
    'words' = words
  )
  if (!is.null(gStats)) {
    ppi = getPPI(org="hs")
    #compute gene-level stats
    p3 = plotGeneStats(gStats, siggs, grps, statName = gStat_name, topN = 5)
    ppi_grps = vissE:::computeMsigGroupPPI(ppi, siggs, grps, gStats) %>%
      activate(nodes) %>% tidygraph::select(group=Group, val=Degree, name=label) %>%
      activate(edges) %>% tidygraph::select(from, to, inferred=Inferred) %>%
      to_split(group, split_by = 'nodes') %>% lapply(as.list) %>% setNames(NULL)

    out$genestats = p3$data
    out$ppi_grps = ppi_grps
  }

  out
}
