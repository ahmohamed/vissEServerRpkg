library(tidyverse)
library(biomaRt)
library(GO.db)
library(GSEABase)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

#----set params----
#creation date
creation_date = as.character(Sys.Date())

#ID types for GSEABase
id_types = c('ensembl_gene_id', 'entrezgene_id', 'external_gene_name', 'uniprot_gn_id')

#output dir
outdir = paste0('inst/extdata/species_gsc')
dir.create(outdir, recursive = TRUE)

#get GO
goids = keys(GO.db, 'GOID')
godb = select(GO.db,  goids, columns(GO.db), 'GOID') |> 
  mutate(
    TERM = gsub(' |,|-', '_', str_to_upper(TERM)),
    TERM = paste0('GO', ONTOLOGY, '_', TERM)
  ) |> 
  column_to_rownames('GOID')

#----build DB for each organism----
#select database
ens = useEnsembl('genes')
organisms = listDatasets(ens)

logs = lapply(organisms$dataset, function(org) tryCatch({
  print(org)
  #build gene ontology for organism
  db = useDataset(org, ens)
  dbattrs = intersect(c('go_id', id_types), listAttributes(db)$name)
  goannotdb = getBM(attributes = dbattrs, 
                # filters = 'chromosome_name', 
                # values = keys(db, 'chromosome_name'), 
                filters = 'ensembl_gene_id', #use for queries that timeout
                values = getBM(attributes = 'ensembl_gene_id', mart = db), 
                mart = db)
  goannotdb = goannotdb |> 
    dplyr::filter(go_id %in% goids) |> 
    mutate(across(everything(), as.character)) |> 
    mutate(across(everything(), ~ ifelse(.x == '', NA_character_, .x)))

  #create ID map
  idmap = goannotdb |> dplyr::select(!go_id) |> distinct()

  #create gene-set collection
  goannotdb = goannotdb |>
    dplyr::select(go_id, ensembl_gene_id) |>
    unique() |> 
    drop_na()
  goannotdb = split(goannotdb, goannotdb$go_id)
  goannotdb = goannotdb[sapply(goannotdb, nrow) >= 5]

  #create gene-sets
  gsc_all = lapply(goannotdb, function (df) {
    go_id = df$go_id[1]
    GeneSet(
      unique(df[, 2]),
      geneIdType = ENSEMBLIdentifier(),
      setName = godb[go_id, 'TERM'],
      shortDescription = godb[go_id, 'DEFINITION'],
      collectionType = BroadCollection('c5', paste0('GO:', godb[go_id, 'ONTOLOGY'])),
      creationDate = creation_date
    )
  })

  #create a GeneSetCollection
  gsc_all = gsc_all[!sapply(gsc_all, is.null)]
  gsc_all = GeneSetCollection(gsc_all)

  #save gscs and ID map
  outfile = file.path(outdir, paste0(gsub('_.+', '_', org), 'gsc', '.rds'))
  saveRDS(gsc_all, outfile)
  outfile = file.path(outdir, paste0(gsub('_.+', '_', org), 'idmap', '.rds'))
  saveRDS(idmap, outfile)
}, error = function(e) e))

# write species information file
organisms |>
  dplyr::rename(Species = "dataset", Description = "description", Version = "version") |>
  mutate(Species = gsub("_.*", "", Species)) |>
  saveRDS("inst/extdata/species_info.rds")

#----create ID maps for human and mouse----
org_cols = c("ENSEMBL", "ENTREZID", "SYMBOL", "UNIPROT")
hsmap = select(org.Hs.eg.db, keys(org.Hs.eg.db, 'SYMBOL'), org_cols, 'SYMBOL')
hsmap = hsmap[, org_cols]
mmmap = select(org.Mm.eg.db, keys(org.Mm.eg.db, 'SYMBOL'), org_cols, 'SYMBOL')
mmmap = mmmap[, org_cols]
colnames(hsmap) = colnames(mmmap) = id_types

#save
saveRDS(hsmap, file.path(outdir, "hsapiens_idmap.rds"))
saveRDS(mmmap, file.path(outdir, "mmusculus_idmap.rds"))

#----ID presence----
id_files = list.files(outdir, recursive = TRUE, full.names = TRUE, pattern = 'idmap.rds')
names(id_files) = gsub('_idmap.rds', '', basename(id_files))
id_present = id_files |> 
  lapply(\(x) {
    data.frame(ID = id_types, Present = id_types %in% colnames(readRDS(x)))
  }) |> 
  bind_rows(.id = 'Species') |> 
  pivot_wider(id_cols = Species, names_from = ID, values_from = Present) |> 
  column_to_rownames('Species') |> 
  as.matrix()
  
# add human and mouse
id_present = rbind(id_present, matrix(TRUE, 2, 4, dimnames = list(c("hsapiens", "mmusculus"), colnames(id_present))))
saveRDS(id_present, file.path("inst/extdata/species_geneid_present.rds"))

#----visualise data availability----
#ID type presence
set.seed(3600)
id_present |> 
  as.data.frame() |> 
  rownames_to_column('Species') |> 
  pivot_longer(!Species, names_to = 'IDType', values_to = 'Present') |> 
  mutate(Present = ifelse(Present, 'Available', 'Unavailable')) |> 
  mutate(
    Species = factor(Species, rownames(id_present)[hclust(dist(id_present + 1))$order]),
    IDType = factor(IDType, colnames(id_present)[hclust(dist(t(id_present) + 1))$order]),
    IDType = case_when(
      IDType == 'external_gene_name' ~ 'Gene symbol',
      IDType == 'ensembl_gene_id' ~ 'Ensembl ID',
      IDType == 'entrezgene_id' ~ 'Entrez ID',
      IDType == 'uniprot_gn_id' ~ 'Uniprot ID'
    )
  ) |>
  ggplot(aes(Species, IDType, fill = Present)) +
  geom_tile() +
  scale_fill_brewer(palette = 'Set2') +
  labs(y = '', title = 'Gene/protein ID types available for each species') +
  vissE::bhuvad_theme() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = 'bottom'
  )

#gene-set size distribution
gsc_files = list.files(outdir, recursive = TRUE, full.names = TRUE, pattern = 'gsc.rds')
tibble(File = gsc_files) |> 
  mutate(Species = gsub('_gsc.rds', '', basename(File))) |> 
  mutate(Length = lapply(File, \(x) sapply(geneIds(readRDS(x)), length))) |> 
  dplyr::select(!File) |> 
  unnest(Length) |> 
  ggplot(aes(group = Species)) +
  geom_density(aes(log2(Length)), colour = alpha('steelblue', 0.1)) +
  labs(x = 'log2(Gene-set length)', y = 'Density', title = 'Gene-set size density across 211 species') +
  vissE::bhuvad_theme()
