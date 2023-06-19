library(ape)
library(ggtree)
library(vegan)
library(TreeDist)
library(phyloseq)
library(tidyverse)
library(glue)
library(qiime2R)
library(ggpubr)

metadata <- read.csv("./ds_metadata.tsv", sep = "\t") %>%
  filter(!row_number() %in% c(1))
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

get_artifact_data <- function(path, ids, extension, metric_list) {
  # Generic import function for artifact data
  artifacts <- list()
  for (id in names(ids)) {
    if (missing(metric_list)) {
      a_path <- glue("{path}/{id}-{extension}.qza")
      artifacts[[id]] <- read_qza(a_path)$data
    } else {
      artifacts[[id]] <- list()
      for (metric in names(metric_list)) {
        a_path <- glue("{path}/{id}/{id}-{extension}{metric_list[[metric]]}")
        artifacts[[id]][[metric]] <- read_qza(glue("{a_path}.qza"))$data
      }
    }
  }
  return(artifacts)
}

known_taxon <- function(row, taxonomy) {
  # Collapse taxonomy into last known taxon
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

merge_with_id <- function(otu_table, taxonomy) {
  # Merge an otu table with a taxonomy table, keeping only identified taxa
  known <- lapply(1:nrow(taxonomy), known_taxon, taxonomy = taxonomy) %>%
    unlist() %>%
    data.frame(row.names = rownames(taxonomy), taxon = .) %>%
    merge(., otu_table, by = 0)
  return(subset(known, select = -c(Row.names)))
}

to_genus_csv <- function(otu_table, taxonomy) {
  # Export a new biom table where the row names have been
  #   replaced with genus-level species identifications where possible
  known <- lapply(1:nrow(taxonomy), genus_level, taxonomy = taxonomy) %>%
    unlist() %>%
    data.frame(row.names = rownames(taxonomy), taxon = .) %>%
    merge(., otu_table, by = 0) %>%
    drop_na()
  rownames(known) <- known$Row.names
  return(subset(known, select = -c(Row.names)))
}

replace_tips <- function(tree, taxonomy_frame) {
  # Map OTU ids to their taxonomic identifications on the tree tips
  taxonomy_frame$known <- lapply(1:nrow(taxonomy_frame), known_taxon,
    taxonomy = taxonomy_frame
  )
  tree$tip.label <- taxonomy_frame$known[tree$tip.label %in%
    rownames(taxonomy_frame)]
  return(tree)
}

metadata_merge_pcoa <- function(metadata, ordination) {
  return(
    ordination %>%
      as.data.frame() %>%
      inner_join(., metadata, by = join_by(x$Vectors.SampleID == y$sample.id))
  )
}

plot_pcoa <- function(pcoa, color_by) {
  return(
    pcoa %>%
      ggplot(aes(
        x = Vectors.PC1, y = Vectors.PC2,
        color = .data[[color_by]]
      )) +
      geom_point()
  )
}

unique_known <- function(otus, identified, classifier) {
  uniques <- identified %>%
    group_by(taxon) %>%
    summarise() %>%
    dim()
  msg <- "% Unique otus:"
  prop <- round((uniques[1] / dim(otus)[1]) * 100, 2)
  return(glue("{classifier} {msg} {prop}"))
}

filter_dm <- function(dm, pattern) {
  # Remove sites from a distance matrix by pattern
  dm <- dm %>%
    as.matrix() %>%
    as.data.frame() %>%
    select(!(matches(pattern))) %>%
    filter(!(grepl(pattern, rownames(.)))) %>%
    as.dist()
  return(dm)
}

filter_meta <- function(metadata, pattern) {
  # Return metadata entries without pattern
  return(metadata %>% filter(!(grepl(pattern, `sample.id`))))
}

sites_x_func <- function(picrust_tsv2) {
  # Reformats picrust's biom output files into a site x function dataframe,
  #   compatible with vegan
  return(picrust_tsv2 %>%
    as.data.frame() %>%
    t() %>%
    `colnames<-`(subset(., grepl("function", rownames(.)))) %>%
    as.data.frame() %>%
    slice(-1) %>%
    mutate_all(as.numeric))
}
