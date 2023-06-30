library(ape)
library(ggtree)
library(TreeSummarizedExperiment)
library(ANCOMBC)
library(vegan)
library(TreeDist)
library(phyloseq)
library(tidyverse)
library(glue)
library(paletteer)
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
  sh = "A_shannon", pi = "A_pielou_e",
  si = "A_simpson", se = "A_simpson_e",
  fa = "APhylo_faith_pd"
)

collapse_tax <- function(taxonomy, level) {
  ranks <- taxonomy %>%
    str_replace_all(".__", "") %>%
    str_split(";") %>%
    unlist()
  if (level > length(ranks)) {
    return(ranks[length(ranks)])
  } else {
    return(ranks[level])
  }
}

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

# Get function annotations
ko <- list()
for (id in names(id_key)) {
  path <- glue("./results/3-FunctionAnnotation/{id}_PICRUST2")
  ko[[id]] <- read_delim(glue("{path}/KO_metagenome_out/pred_metagenome_unstrat.tsv"))
}

rel_abund <- function(abs_abund, first_col) {
    rel_abund <- data.frame(first_col = abs_abund[1])
    for (col in 2:ncol(abs_abund)) {
        rel_abund[colnames(abs_abund[col])] <- abs_abund[col] / sum(abs_abund[col])
    }
    return(rel_abund)
}

known_taxon <- function(row, taxonomy, level) {
    # Collapse taxonomy into last known taxon or specified taxonomic rank
    if (missing(level)) {
        known_rank <- 7
        while (is.na(taxonomy[row, known_rank]) && known_rank != 1) {
            known_rank <- known_rank - 1
        }
        return(taxonomy[row, known_rank])
    } else {
        return(taxonomy[row, level])
    }
}

combine_freqs <- function(freq_list, sum_by) {
  combined <- bind_rows(freq_list) %>%
    arrange(.[["sum_by"]]) %>%
    group_by(sum_by) %>%
    summarise(across(everything(), sum)) %>%
    mutate_all(~ replace(., is.na(.), 0))
  return(combined)
}

genus_level <- function(row, taxonomy) {
  if (is.na(taxonomy[row, 7])) {
    return(taxonomy[row, 6])
  }
  return(taxonomy[row, 7])
}

merge_with_id <- function(otu_table, taxonomy, level) {
    # Merge an otu table with a taxonomy table, keeping only identified taxa
    known <- lapply(1:nrow(taxonomy), known_taxon,
        taxonomy = taxonomy,
        level = level
    ) %>%
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

replace_tips <- function(tree, taxonomy_frame, level) {
  # Map OTU ids to their taxonomic identifications on the tree tips
  known <- lapply(1:nrow(taxonomy_frame), known_taxon,
    taxonomy = taxonomy_frame,
    level = level
  )
  tree$tip.label <- known[tree$tip.label %in%
    rownames(taxonomy_frame)]
  return(tree)
}

metadata_merge_pcoa <- function(metadata, ordination, functions) {
  if (missing(functions)) {
    return(
      ordination %>%
        as.data.frame() %>%
        inner_join(., metadata, by = join_by(x$Vectors.SampleID == y$sample.id))
    )
  } else {
    return(ordination %>%
      as.data.frame() %>%
      mutate(sample.id = rownames(.)) %>%
      inner_join(metadata))
  }
}

plot_pcoa <- function(pcoa, color_by, functions, title, subtitle = NULL) {
  if (missing(functions)) {
    x <- "Vectors.PC1"
    y <- "Vectors.PC2"
  } else {
    x <- "V1"
    y <- "V2"
  }
  type <- pcoa[["Type"]]
  return(
    pcoa %>%
      ggplot(aes(
        x = .data[[x]], y = .data[[y]],
        color = .data[[color_by]]
      )) +
      geom_point(
        aes(shape = .data[["Type"]]),
        size = 2,
        stroke = 1
      )
      +
      scale_color_paletteer_d("pals::glasbey") +
      labs(x = "PC1", y = "PC2", title = title, subtitle = subtitle)
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

filter_dm <- function(dm, keep) {
  # Keep only specified sites in a distance matrix
  dm <- dm %>%
    as.matrix() %>%
    as.data.frame() %>%
    select(matches(keep)) %>%
    filter(grepl(paste(keep, collapse = "|"), rownames(.))) %>%
    as.dist()
  return(dm)
}

filter_meta <- function(metadata, keep) {
  # Keep only specified sites in metadata with pattern
  return(metadata %>% filter(grepl(paste(keep, collapse = "|"), `sample.id`)))
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

import_ancom <- function(result, path) {
  return(read_csv(glue("{path}/export/{result}")) %>%
    select(.data = ., !`(Intercept)`))
}
