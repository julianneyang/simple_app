
library(tidyverse)

smt_full <- read.delim("results/Bracken/SMT_Maaslin3/Full/Sequencing_Depth_Site_Sex_Genotype/significant_results.tsv")

slc_full <- read.delim("results/Bracken/SPONT_Maaslin3/Full/Sequencing_Depth_Site_Sex_Genotype/all_results.tsv") 

filtered_slc_full_abundance <- slc_full %>% 
  filter(metadata=="Genotype") %>% 
  filter(null_hypothesis!=0)

# Show only DEG or pathway that overlap -

make_overlap_plot <- function(het_input,
                              mut_input,
                              threshold = 1,
                              type = c("DEG", "GSEA"),
                              padj_cutoff = 0.25) {
  
  library(dplyr)
  library(tidyr)
  library(plotly)
  
  type <- match.arg(type)
  # Column mapping
  id_col   <- if (type == "DEG") "external_gene_name" else "pathway"
  coef_col <- if (type == "DEG") "log2FoldChange"     else "NES"
  desc_col <- "description"   # optional; may not exist for some files
  
  # helper to read either a path or accept a pre-loaded dataframe
  read_input <- function(x) {
    if (is.character(x) && file.exists(x)) {
      df <- read.csv(x, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (is.data.frame(x)) {
      df <- x
    } else {
      stop("het_input / mut_input must be a file path or a data.frame")
    }
    # If id_col is not a column, try to recover from rownames
    if (!(id_col %in% names(df))) df[[id_col]] <- rownames(df)
    df
  }
  
  het <- read_input(het_input) %>% mutate(value = "HET")
  mut <- read_input(mut_input) %>% mutate(value = "MUT")
  
  # Are padj / description present?
  has_padj <- "padj" %in% names(het) || "padj" %in% names(mut)
  has_desc <- desc_col %in% names(het) || desc_col %in% names(mut)
  
  # If padj exists, create list of significant features from each file
  if (has_padj) {
    het_sig <- het %>% filter(padj < padj_cutoff) %>% pull(!!sym(id_col))
    mut_sig <- mut %>% filter(padj < padj_cutoff) %>% pull(!!sym(id_col))
    combined_sig <- union(het_sig, mut_sig)
  } else {
    combined_sig <- NULL
  }
  
  # Combine and optionally keep only combined significant features
  full <- bind_rows(het, mut)
  if (!is.null(combined_sig)) {
    full <- full %>% filter(.data[[id_col]] %in% combined_sig)
  }
  
  # Build HET and MUT side-by-side by selecting & renaming
  # Use mutate() to create HET/MUT columns from coef_col; also keep padj & description if present
  het_df <- full %>%
    filter(value == "HET") %>%
    mutate(HET = .data[[coef_col]],
           padj_HET = if (has_padj) .data$padj else NA_real_,
           description_HET = if (has_desc) .data[[desc_col]] else NA_character_) %>%
    select(all_of(c(id_col, "HET",
                    if (has_padj) "padj_HET" else NULL,
                    if (has_desc) "description_HET" else NULL)))
  
  mut_df <- full %>%
    filter(value == "MUT") %>%
    mutate(MUT = .data[[coef_col]],
           padj_MUT = if (has_padj) .data$padj else NA_real_,
           description_MUT = if (has_desc) .data[[desc_col]] else NA_character_) %>%
    select(all_of(c(id_col, "MUT",
                    if (has_padj) "padj_MUT" else NULL,
                    if (has_desc) "description_MUT" else NULL)))
  
  # Inner join ensures only features with both HET & MUT remain
  full_wide <- inner_join(het_df, mut_df, by = id_col)
  
  # Keep only same-direction features and those meeting the absolute threshold
  full_wide <- full_wide %>%
    filter(!is.na(HET) & !is.na(MUT)) %>%
    filter((HET > 0 & MUT > 0) | (HET < 0 & MUT < 0)) %>%
    filter((abs(HET) >= threshold) | (abs(MUT) >= threshold))
  
  if (nrow(full_wide) == 0) {
    # return an empty informative plotly object
    return(plot_ly() %>% layout(title = "No features passed filters"))
  }
  
  # Combine description if available
  if (has_desc) {
    # description from HET preferred, else from MUT
    if ("description_HET" %in% names(full_wide) || "description_MUT" %in% names(full_wide)) {
      full_wide <- full_wide %>%
        mutate(description = coalesce(description_HET, description_MUT)) %>%
        select(-any_of(c("description_HET", "description_MUT")))
    } else {
      full_wide$description <- NA_character_
    }
  } else {
    full_wide$description <- NA_character_
  }
  
  # compute ordering by mean coefficient
  full_wide <- full_wide %>%
    mutate(mean_coef = (HET + MUT) / 2) %>%
    arrange(mean_coef)
  
  ordered_features <- full_wide[[id_col]]
  
  # long format suitable for plotly bars
  full_long <- full_wide %>%
    select(all_of(c(id_col, "HET", "MUT", "description",
                    if (has_padj) c("padj_HET","padj_MUT") else NULL))) %>%
    pivot_longer(cols = c("HET", "MUT"), names_to = "value", values_to = "coef") %>%
    mutate(feature = factor(.data[[id_col]], levels = ordered_features))
  
  # Build plotly grouped horizontal bar chart
  hover_text <- if (has_padj) {
    paste0("<b>", full_long$feature, "</b><br>",
           "Condition: ", full_long$value, "<br>",
           coef_col, ": ", round(full_long$coef, 2), "<br>",
           "Description: ", full_long$description, "<br>",
           "padj (if present): ", signif(ifelse(full_long$value == "HET", full_long$padj_HET, full_long$padj_MUT), 3))
  } else {
    paste0("<b>", full_long$feature, "</b><br>",
           "Condition: ", full_long$value, "<br>",
           coef_col, ": ", round(full_long$coef, 2), "<br>",
           "Description: ", full_long$description)
  }
  
  fig <- plot_ly(
    data = full_long,
    x = ~coef,
    y = ~feature,
    color = ~value,
    colors = c("steelblue", "tomato"),
    type = "bar",
    orientation = "h",
    text = hover_text,
    hoverinfo = "text"
  ) %>%
    layout(
      barmode = "group",
      yaxis = list(title = "", automargin = TRUE),
      xaxis = list(title = ifelse(type == "DEG", "log2 Fold Change", "Normalized Enrichment Score (NES)")),
      hoverlabel = list(bgcolor = "white")
    )
  
  # return both objects
  return(list(
    plot = fig,
    table = full_wide
  ))
}