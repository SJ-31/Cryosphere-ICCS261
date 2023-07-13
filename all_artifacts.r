library(microseq)
library(caret)
library(randomForest)
library(ape)
library(ggridges)
library(ggpattern)
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
  filter(!row_number() %in% c(1)) # Read in the metadata

id_key <- list(
  # Prepare the site list
  BhM = "Bihor mountains", BrH = "Barrow mountain high",
  BrL = "Barrow mountain low", CaS = "Catriona snow", CeY = "Central Yakutia",
  GrI = "Greenland ice", EaI = "East Iceland glaciers", CrC = "Cryoconite",
  HaS = "Hailstone", NzS = "New Zealand soil", SvG = "Sverdrup glacier",
  StR = "Storglaciaren", ViS = "Villum station"
)
beta_metrics <- list(
  # List of beta metrics
  bc = "B_braycurtis", ja = "B_jaccard",
  uu = "BPhylo_unweighted_unifrac",
  wn = "BPhylo_weighted_normalized_unifrac",
  wu = "BPhylo_weighted_unifrac"
)
alpha_metrics <- list(
  # list of alpha metrics
  sh = "A_shannon", pi = "A_pielou_e",
  si = "A_simpson", se = "A_simpson_e",
  fa = "APhylo_faith_pd"
)

collapse_tax <- function(taxonomy, level) {
  # Collapses taxonomy dataframe into a specified taxonomic rank
  ranks <- taxonomy %>%
    str_replace_all(".__", "") %>%
    str_split(";") %>%
    unlist() # The R list structure is actually like a dictionary so this funciton flattens it
  if (level > length(ranks)) {
    return(ranks[length(ranks)])
  } else {
    return(ranks[level])
  }
}

count_identified <- function(taxonomy, name) {
  # Count the number of identified taxa
  total_ranks <- taxonomy %>%
    lapply(., dim) %>%
    unlist(use.names = FALSE) %>%
    sum()
  total_ranks <- total_ranks - (length(taxonomy) * 7)
  identified <- taxonomy %>%
    lapply(., is.na) %>%
    lapply(., which) %>%
    unlist() %>%
    length()
  ratio <- 1 - (identified / (total_ranks * 7))
  glue("{name}: {ratio}")
}

get_artifact_data <- function(path, ids, extension, metric_list) {
  # Generic import function for artifact data
  artifacts <- list()
  for (id in names(ids)) {
    if (missing(metric_list)) {
      # glue is the equivalent of an F-string
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
  ko[[id]] <- read_delim(glue("{path}/pathways_out/path_abun_unstrat.tsv"))
}

rel_abund <- function(abs_abund, first_col) {
  # Calculate relative abundance from absolute abundance
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
  # Combine frequency tables from different sites into one, summing up the frequencies
  combined <- bind_rows(freq_list) %>%
    arrange(.[["sum_by"]]) %>%
    group_by(sum_by) %>%
    summarise(across(everything(), sum)) %>%
    mutate_all(~ replace(., is.na(.), 0))
  # replace na with 0
  return(combined)
}

genus_level <- function(row, taxonomy) {
  # Return genus-level identifications
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
  # Merge metadata with pcoa results
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

plot_pcoa <- function(pcoa, color_by, functions) {
  # Plot pcoa graph
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
        size = 2, # Specifications for points
        stroke = 1
      )
      +
      scale_color_paletteer_d("pals::glasbey") +
      theme(
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)
      )
  )
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
  #   compatible with vegan's vegdist function
  return(picrust_tsv2 %>%
    as.data.frame() %>%
    t() %>%
    `colnames<-`(subset(., grepl("function", rownames(.)))) %>%
    as.data.frame() %>%
    dplyr::slice(-1) %>%
    mutate_all(as.numeric))
}

import_ancom <- function(result, path) {
  return(read_csv(glue("{path}/export/{result}")) %>%
    select(.data = ., !`(Intercept)`))
}

ancombc_select <- function(ancombc_results, result, tax_level, unwanted) {
  # Select results type from ancombc results object in long format
  select <- ancombc_results %>%
    select(c(1, grep(result, colnames(ancombc_results))))
  if (!(is.na(tax_level)) || !(is.na(unwanted))) { # Select the specified taxonomic level and ignore unwanted
    select <- select %>%
      filter(grepl(tax_level, .data$taxon)) %>%
      filter(!(grepl(paste(unwanted, collapse = "|"), .data$taxon)))
  } else {
    tax_level <- NaN
  }
  selected <- select %>%
    mutate(taxon = str_replace(taxon, glue("{tax_level}:"), "")) %>%
    `colnames<-`(str_replace(colnames(.), result, "")) %>%
    pivot_longer(., -taxon)
  colnames(selected)[which(names(selected) == "value")] <- result
  return(selected)
}

prepare_abc_lfc <- function(abc_results, var, results_type, rank, wrong_tax) {
  # Prepare ancombc log fold change results for plotting
  old <- glue("lfc_{var} se_{var} diff_{var}") %>%
    strsplit(" ") %>%
    unlist()
  new <- c("lfc", "se", "diff")
  abc <- abc_results[[results_type]] %>%
    select(-(grep("Intercept", colnames(abc_results[[results_type]]))))
  diff_abund <- ancombc_select(abc, glue("diff_{var}"), rank, wrong_tax)
  se <- ancombc_select(abc, glue("se_{var}"), rank, wrong_tax)
  lfc <- ancombc_select(abc, glue("lfc_{var}"), rank, wrong_tax) %>%
    merge(se, by = c("name", "taxon")) %>%
    merge(diff_abund, by = c("name", "taxon")) %>%
    filter((!!as.symbol(glue("diff_{var}"))) == TRUE) %>%
    rename_with(~new, all_of(old))
  return(lfc)
  # Don't want to show taxa that don't have statistically
  # signifcant differences in log fold change
}

quantile_filter <- function(lfc_table, cutoff) {
  # Filter log fold change so that sites above the 3rd and below the 1st quartile remain
  upper <- lfc_table$lfc %>% quantile(1 - cutoff)
  lower <- lfc_table$lfc %>% quantile(cutoff)
  highs_lows <- lfc_table[lfc_table$lfc > upper | lfc_table$lfc < lower, ]
  return(highs_lows)
}


sum_by_site <- function(freq_table, id_key, id_col, unwanted) {
  # Add up frequencies for different samples of the same site
  summed <- sapply(names(id_key), function(x) {
    rowSums(freq_table[, grep(x, colnames(freq_table)), drop = FALSE])
  }) %>%
    as_tibble() %>%
    mutate(identifier = freq_table[[id_col]]) %>%
    relocate(identifier) %>%
    filter(!(grepl("[0-9]", identifier))) %>%
    # Remove uncharacterized taxa
    filter(!(grepl(paste(unwanted, collapse = "|"), identifier))) %>%
    rel_abund(., identifier) %>%
    pivot_longer(., -identifier) %>%
    mutate(name = id_key[.data$name] %>% unlist(use.names = FALSE))
  return(summed)
}

abc_lfc_plot <- function(abc_lfc) {
  # Plot ancombc log fold change
  plot <- abc_lfc %>%
    ggplot(aes(x = name, y = lfc, fill = taxon)) +
    geom_bar(
      stat = "identity",
      position = position_dodge()
    ) +
    geom_errorbar(
      aes(
        ymin = lfc - se,
        ymax = lfc + se
      ),
      width = .2,
      position = position_dodge(.9)
    ) +
    labs(x = "Site", y = "Log fold change")
  return(plot)
}
