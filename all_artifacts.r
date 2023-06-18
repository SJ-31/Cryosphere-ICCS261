library(ape)
library(TreeDist)
library(phyloseq)
library(tidyverse)
library(glue)
library(qiime2R)

id_key <- list(
  BhM = "Bihor mountains", BrH = "Barrow mountain high",
  BrL = "Barrow mountain low", CaS = "Catriona snow", CeY = "Central Yakutia",
  GrI = "Greenland ice", EaI = "East Iceland glaciers", CrC = "Cryoconite",
  HaS = "Hailstone", NzS = "New Zealand soil", SvG = "Sverdrup glacier",
  StR = "Storglaciaren", ViS = "Villum station"
)
beta_metrics <- list(
  bc = "B_braycurtis", ja = "B_jaccard",
  uu = "BPhylo_unweighted_unifrac",
  wn = "BPhylo_weighted_normalized_unifrac",
  wu = "BPhylo_weighted_unifrac"
)
alpha_metrics <- list(
  sh = "shannon", pi = "pielou_e",
  si = "simpson", se = "simpson_e"
)

get_artifact_data <- function(path, ids, extension) {
  # Generic import function for artifact data
  artifacts <- list()
  for (id in names(ids)) {
    a_path <- glue("{path}/{id}-{extension}.qza")
    artifacts[[id]] <- read_qza(a_path)$data
  }
  return(artifacts)
}

known_taxon <- function(row, taxonomy) {
  known_rank <- 7
  while (is.na(taxonomy[row, known_rank]) && known_rank != 1) {
    known_rank <- known_rank - 1
  }
  return(taxonomy[row, known_rank])
}

genus_level <- function(row, taxonomy) {
  if (is.na(taxonomy[row, 7])) {
    return(taxonomy[row, 6])
  }
  return(taxonomy[row, 7])
}

to_genus_csv <- function(otu_table, taxonomy, file_name) {
  # Export a new biom table where the row names have been
  #   replaced with genus-level species identifications where possible
  known <- lapply(1:nrow(taxonomy), genus_level, taxonomy = taxonomy) %>%
    unlist() %>%
    data.frame(rownames = rownames(taxonomy), taxon = .)
  write.csv(known %>% drop_na(), file = file_name)
}

replace_tips <- function(tree, taxonomy_frame) {
  # Map OTU ids to their taxonomic identifications on the tree tips
  taxonomy_frame$known <- lapply(seq_le(nrow(taxonomy_frame)), known_taxon,
    taxonomy = taxonomy_frame
  )
  tree$tip.label <- taxonomy_frame$known[tree$tip.label %in%
    rownames(taxonomy_frame)]
  return(tree)
}