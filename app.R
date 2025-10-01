library(shiny)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)
library(here)
library(DT)
library(plotly)
library(stringi)
library(nlme)


# ---- Load data ----

# ---- Global functions ----
plot_histology <- function(csv_filepath, title_string, subtitle_string, stat_comparisons,
                           color_vector = c("WT"="black","HET"= "navy","MUT"="firebrick")) {
  histo <- readr::read_csv(csv_filepath)
  histo$Genotype <- factor(histo$Genotype, levels=c("WT", "HET", "MUT"))
  
  ggplot(data=histo,aes(x=Genotype,y=Score, color=Genotype)) + 
    stat_summary(aes(x=Genotype, y=Score), fun=median, geom="crossbar", colour="black")+
    geom_beeswarm(cex = 3,priority = "density",size=3)+
    scale_color_manual(values=color_vector)+
    theme_cowplot(12) +
    theme(legend.position = "none")+
    ggtitle({{title_string}})+
    labs(subtitle = {{subtitle_string}})+
    stat_compare_means(comparisons = {{stat_comparisons}},method = "wilcox")+
    ylab("Score")+
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust=0.5))
  
}

plot_fitc <- function(csv_filepath, title_string, subtitle_string, stat_comparisons,
                      color_vector = c("WT"="black","HET"= "navy","MUT"="firebrick")) {
  df <- readr::read_csv(csv_filepath)
  df$Genotype <- factor(df$Genotype, levels=c("WT", "HET", "MUT"))
  
  ggplot(df, aes(x = Genotype, y = Plasma_FITC, fill = Genotype)) +
    geom_boxplot(alpha=0.5) +
    scale_fill_manual(values=color_vector) +
    geom_jitter(width = 0.2, alpha=0.8) +
    theme_cowplot(12) +
    theme(legend.position = "none")+
    ggtitle({{title_string}})+
    labs(subtitle = {{subtitle_string}})+
    stat_compare_means(comparisons = {{stat_comparisons}},method = "t.test")+
    ylab("Score")+
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust=0.5))
  
}

plot_avg_trajectory <- function(df, subtitle_string, title_string,                       
                                color_vector = c("WT"="black","HET"= "navy","MUT"="firebrick")) {
  
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  
  # Reshape to long format
  df_long <- df %>%
    pivot_longer(
      cols = starts_with("X"), 
      names_to = "Time", 
      values_to = "Value"
    ) %>%
    mutate(
      Time = as.numeric(gsub("X|_hours", "", Time)) # convert to numeric hours
    )
  
  summary(df_long$Value)
  df_long <- drop_na(df_long)
  
  #df_long <- replace_na(df_long, list(Value = 9))
  df_long$Genotype <- factor(df_long$Genotype, levels=c("WT", "HET","MUT"))
  
  lm <- lme(Value ~ Sex + Time*Genotype, data = df_long, random = ~ 1| MouseID)
  print(summary(lm))
  
  # Summarize by Genotype and Time
  df_summary <- df_long %>%
    group_by(Genotype, Time) %>%
    summarise(
      mean_value = mean(Value, na.rm = TRUE),
      se = sd(Value, na.rm = TRUE) / sqrt(n())
    )
  
  # Plot average ± SE
  ggplot(df_summary, aes(x = Time, y = mean_value, color = Genotype, group = Genotype)) +
    geom_line(size = 1.2) +
    scale_color_manual(values={{color_vector}})+
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean_value - se, ymax = mean_value + se), width = 2) +
    labs(x = "Time (hours)", y = "Mean Value ± SE", title = "Average Trajectory by Genotype") +
    theme_cowplot(12) +
    ggtitle({{title_string}})+
    labs(subtitle = {{subtitle_string}})+
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust=0.5))
}

# Helper function: read, preprocess, and make GSEA Manhattan plot
set.seed(123)

make_gsea_plot <- function(path) {
  df <- read.csv(here(path)) %>%
    as.data.frame() %>%
    mutate(logFDR = -log10(padj))
  
  plot_ly(
    data = df,
    x = ~NES,
    y = ~logFDR,
    text = ~paste0(
      "Pathway: ", pathway, "<br>",
      "NES: ", round(NES, 2), "<br>",
      "FDR: ", signif(padj, 3)
    ),
    type = "scatter",
    mode = "markers",
    marker = list(
      size = 10,
      color = ~logFDR,
      colorscale = "Viridis",
      showscale = TRUE
    )
  ) %>%
    layout(
      xaxis = list(title = "Normalized Enrichment Score (NES)"),
      yaxis = list(title = "-log10(padj)"),
      hovermode = "closest"
    )
}

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




make_comparison_plot_smt <- function(
    gsea_file_smt,
    gsea_file_spont,
    threshold = 1,
    filter_concordant = FALSE,
    apply_threshold = FALSE
) {
  # 1. Read and format both datasets
  deg_smt <- read.csv(gsea_file_smt, row.names = 1) %>%
    mutate(value = "SMT", mean_coef = NES) %>%
    dplyr::select(pathway, mean_coef, value, padj)
  
  deg_spont <- read.csv(gsea_file_spont) %>%
    mutate(value = "SPONT", padj = padj_MUT) %>%
    dplyr::select(pathway, mean_coef, value, padj)
  
  # 2. Combine them
  full <- rbind(deg_smt, deg_spont)
  
  # 3. Keep significant features in either dataset
  significant_smt <- deg_smt %>% filter(padj < 0.25) %>% pull(pathway)
  significant_spont <- deg_spont %>% pull(pathway)
  combined_significant_features <- union(significant_smt, significant_spont)
  full <- full %>% filter(pathway %in% combined_significant_features)
  
  # 4. Require both conditions to exist
  full_filtered <- full %>%
    group_by(pathway) %>%
    filter(all(c("SMT", "SPONT") %in% value)) %>%
    ungroup()
  
  # 5. Optionally filter for concordance
  if (filter_concordant) {
    full_filtered <- full_filtered %>%
      group_by(pathway) %>%
      filter({
        smt_coef <- mean_coef[value == "SMT"]
        spont_coef <- mean_coef[value == "SPONT"]
        length(smt_coef) == 1 && length(spont_coef) == 1 &&
          ((smt_coef > 0 & spont_coef > 0) | (smt_coef < 0 & spont_coef < 0))
      }) %>%
      ungroup()
  }
  
  # 6. Optionally filter by effect size threshold
  if (apply_threshold) {
    full_filtered <- full_filtered %>%
      group_by(pathway) %>%
      filter({
        smt_coef <- mean_coef[value == "SMT"]
        spont_coef <- mean_coef[value == "SPONT"]
        abs(smt_coef) >= threshold | abs(spont_coef) >= threshold
      }) %>%
      ungroup()
  }
  
  # 7. Calculate mean effect size per pathway for ordering
  ordered_pathways <- full_filtered %>%
    group_by(pathway) %>%
    summarise(avg_coef = mean(mean_coef, na.rm = TRUE)) %>%
    arrange(avg_coef) %>%
    pull(pathway)
  
  # 8. Prepare data for plotting
  full_plot <- full_filtered %>%
    mutate(pathway = factor(pathway, levels = ordered_pathways))
  
  # 9. Create the interactive plot
  fig <- plot_ly(
    data = full_plot,
    x = ~mean_coef,
    y = ~pathway,
    color = ~value,
    colors = c("steelblue", "tomato"),
    type = "bar",
    orientation = "h",
    hoverinfo = "text",
    text = ~paste0(
      "<b>", pathway, "</b><br>",
      "Condition: ", value, "<br>",
      "Effect size: ", round(mean_coef, 2), "<br>",
      "padj: ", signif(padj, 3)
    )
  ) %>%
    layout(
      barmode = "group",
      yaxis = list(title = "", categoryorder = "array", categoryarray = ordered_pathways),
      xaxis = list(title = "Effect size"),
      hoverlabel = list(bgcolor = "white")
    )
  
  # 10. Wide-format table
  full_wide <- full_filtered %>%
    dplyr::select(pathway, value, mean_coef, padj) %>%
    tidyr::pivot_wider(names_from = value, values_from = c(mean_coef, padj))
  print(head(full_wide))
  
  return(list(plot = fig, table = full_wide))
}

make_comparison_plot <- function(
    gsea_file_1,
    gsea_file_spont,
    threshold = 1,
    filter_concordant = FALSE,
    apply_threshold = FALSE
) {
  # 1. Read and format both datasets
  deg_smt <- data.frame(gsea_file_1) %>%
    mutate(value = "CURRENT") %>%
    dplyr::select(pathway, mean_coef, value, padj_MUT,padj_HET)
  
  deg_spont <- read.csv(gsea_file_spont) %>%
    mutate(value = "SPONT") %>%
    dplyr::select(pathway, mean_coef, value, padj_MUT,padj_HET)
  
  # 2. Combine them
  full <- rbind(deg_smt, deg_spont)
  
  # 3. Keep significant features in either dataset
  significant_smt <- deg_smt %>% pull(pathway)
  significant_spont <- deg_spont %>% pull(pathway)
  combined_significant_features <- union(significant_smt, significant_spont)
  full <- full %>% filter(pathway %in% combined_significant_features)
  
  # 4. Require both conditions to exist
  full_filtered <- full %>%
    group_by(pathway) %>%
    filter(all(c("CURRENT", "SPONT") %in% value)) %>%
    ungroup()
  
  # 5. Optionally filter for concordance
  if (filter_concordant) {
    full_filtered <- full_filtered %>%
      group_by(pathway) %>%
      filter({
        smt_coef <- mean_coef[value == "CURRENT"]
        spont_coef <- mean_coef[value == "SPONT"]
        length(smt_coef) == 1 && length(spont_coef) == 1 &&
          ((smt_coef > 0 & spont_coef > 0) | (smt_coef < 0 & spont_coef < 0))
      }) %>%
      ungroup()
  }
  
  # 6. Optionally filter by effect size threshold
  if (apply_threshold) {
    full_filtered <- full_filtered %>%
      group_by(pathway) %>%
      filter({
        smt_coef <- mean_coef[value == "CURRENT"]
        spont_coef <- mean_coef[value == "SPONT"]
        abs(smt_coef) >= threshold | abs(spont_coef) >= threshold
      }) %>%
      ungroup()
  }
  
  # 7. Calculate mean effect size per pathway for ordering
  ordered_pathways <- full_filtered %>%
    group_by(pathway) %>%
    summarise(avg_coef = mean(mean_coef, na.rm = TRUE)) %>%
    arrange(avg_coef) %>%
    pull(pathway)
  
  # 8. Prepare data for plotting
  full_plot <- full_filtered %>%
    mutate(pathway = factor(pathway, levels = ordered_pathways))
  
  # 9. Create the interactive plot
  fig <- plot_ly(
    data = full_plot,
    x = ~mean_coef,
    y = ~pathway,
    color = ~value,
    colors = c("steelblue", "tomato"),
    type = "bar",
    orientation = "h",
    hoverinfo = "text",
    text = ~paste0(
      "<b>", pathway, "</b><br>",
      "Condition: ", value, "<br>",
      "Effect size: ", round(mean_coef, 2), "<br>",
      "padj_MUT: ", signif(padj_MUT, 3),
      "padj_HET: ", signif(padj_HET, 3)
    )
  ) %>%
    layout(
      barmode = "group",
      yaxis = list(title = "", categoryorder = "array", categoryarray = ordered_pathways),
      xaxis = list(title = "Effect size"),
      hoverlabel = list(bgcolor = "white")
    )
  
  # 10. Wide-format table
  full_wide <- full_filtered %>%
    dplyr::select(pathway, value, mean_coef, padj_HET,padj_MUT) %>%
    tidyr::pivot_wider(names_from = value, values_from = c(mean_coef, padj_HET,padj_MUT))
  print(head(full_wide))
  
  return(list(plot = fig, table = full_wide))
}

# ---- UI ----
ui <- fluidPage(
  titlePanel("Aggregated Data from the SLC SI Project"),
  
  tabsetPanel(id = "tabs",
              tabPanel("SLC_Spontaneous",
                       fluidRow(
                         h2("Phenotype Results"),
                         column(4, plotOutput("plot_slc_spont_1")),
                         column(4, plotOutput("plot_slc_spont_2")),
                         column(4, plotOutput("plot_slc_spont_3")),
                         h2("RNA-sequencing Results"),
                         br(),
                         h3("DESeq2: MUT vs WT"),
                         textOutput("filepath01"),
                         DTOutput("preview01"),
                         h3("DESeq2: HET vs WT"),
                         textOutput("filepath02"),
                         DTOutput("preview02"),
                         br(),
                         h3("GSEA: MUT vs WT"),
                         textOutput("filepath03"),
                         DTOutput("preview03"),
                         plotlyOutput("gsea_manhattan_mut", height = "600px"),
                         h3("GSEA: HET vs WT"),
                         textOutput("filepath05"),
                         DTOutput("preview05"),
                         plotlyOutput("gsea_manhattan_het", height = "600px"),
                         br(),
                         h3("DEG with concordant directionality between HET and MUT, and padj < 0.25 in at least one"),
                         sliderInput("threshold", "Absolute log2FC threshold:", 
                                     min = 0, max = 3, value = 1, step = 0.1),
                         DTOutput("overlap_table"),
                         downloadButton("download_full_wide", "Download Concordant DEG CSV"),
                         plotlyOutput("plot_slc_spont_4"),
                         h3("Pathways with concordant directionality between HET and MUT, and padj < 0.25 in at least one"),
                         sliderInput("gsea_threshold", "Absolute log2FC threshold:", 
                                     min = 0, max = 3, value = 1, step = 0.1),
                         DTOutput("path_overlap_table"),
                         downloadButton("download_path_wide", "Download Concordant Pathways CSV"),
                         plotlyOutput("plot_slc_spont_5")
                       )
                       
              ),
              
              tabPanel("SMT_Negative",
                       fluidRow(
                         h2("Phenotype Results"),
                         column(6, plotOutput("plot_smt_negative_1")),
                         column(6, plotOutput("plot_smt_negative_2")),
                         br(),
                         h3("SMT Shotgun Data"),
                         column(6, plotOutput("plot_smt_negative_3")),
                         column(6, plotOutput("plot_smt_negative_4")),
                         br(),
                         h3("DESeq2: MUT vs WT"),
                         textOutput("filepath07"),
                         DTOutput("preview07"),
                         br(),
                         h3("GSEA: MUT vs WT"),
                         textOutput("filepath08"),
                         DTOutput("preview08"),
                         plotlyOutput("SMT_Neg_gsea_manhattan", height = "600px"),
                         h3("Concordance with SPONT FITC"),
                         checkboxInput("filter_concordant", "Filter for concordant direction", value = FALSE),
                         checkboxInput("apply_threshold", "Filter by effect size", value = FALSE),
                         sliderInput("threshold", "Effect size threshold", min = 0, max = 3, step = 0.1, value = 1),
                         plotlyOutput("comparison_plot"),
                         DTOutput("comparison_table"),
                         downloadButton("download_comparison")
                       )
              ),
              
              tabPanel("SLC_TL1A",
                       fluidRow(
                         column(6, plotOutput("plot_tl1a_1")),
                         column(6, plotOutput("plot_tl1a_2")),
                         h2("RNA-sequencing Results"),
                         br(),
                         h3("DESeq2: MUT vs WT"),
                         textOutput("filepath13"),
                         DTOutput("preview13"),
                         h3("DESeq2: HET vs WT"),
                         textOutput("filepath12"),
                         DTOutput("preview12"),
                         h3("GSEA: MUT vs WT"),
                         textOutput("filepath14"),
                         DTOutput("preview14"),
                         h3("GSEA: HET vs WT"),
                         textOutput("filepath15"),
                         DTOutput("preview15"),
                         h3("DEG with concordant directionality between HET and MUT, and padj < 0.25 in at least one"),
                         sliderInput("threshold", "Absolute log2FC threshold:", 
                                     min = 0, max = 3, value = 1, step = 0.1),
                         DTOutput("tl1a_overlap_table"),
                         plotlyOutput("plot_tl1a_4"),
                         h3("Pathways with concordant directionality between HET and MUT, and padj < 0.25 in at least one"),
                         sliderInput("gsea_threshold", "Absolute log2FC threshold:", 
                                     min = 0, max = 3, value = 1, step = 0.1),
                         DTOutput("tl1a_path_overlap_table"),
                         plotlyOutput("plot_tl1a_5"),
                         h3("Concordance with SPONT FITC"),
                         checkboxInput("tl1a_filter_concordant", "Filter for concordant direction", value = FALSE),
                         checkboxInput("tl1a_apply_threshold", "Filter by effect size", value = FALSE),
                         sliderInput("tl1a_threshold", "Effect size threshold", min = 0, max = 3, step = 0.1, value = 1),
                         DTOutput("tl1a_comparison_table"),
                         downloadButton("tl1a_download_comparison"),
                         plotlyOutput("tl1a_comparison_plot")
                       )
              ),
              
              tabPanel("SLC_Indomethacin",
                       fluidRow(
                         column(6, plotOutput("plot_slc_indo_1")),
                         column(6, plotOutput("plot_slc_indo_2"))
                       )
              ),
              
              tabPanel("SLC_HFD",
                       fluidRow(
                         h2("Phenotype Results"),
                         column(6, plotOutput("plot_slc_hfd_1")),
                         column(6, plotOutput("plot_slc_hfd_2")),
                         h2("RNA-sequencing Results"),
                         br(),
                         h3("DESeq2: MUT vs WT"),
                         textOutput("filepath10"),
                         DTOutput("preview10"),
                         h3("DESeq2: HET vs WT"),
                         textOutput("filepath16"),
                         DTOutput("preview16"),
                         h3("GSEA: MUT vs WT"),
                         textOutput("filepath11"),
                         DTOutput("preview11"),
                         h3("GSEA: HET vs WT"),
                         textOutput("filepath17"),
                         DTOutput("preview17"),
                         br(),
                         h3("DEG with concordant directionality between HET and MUT, and padj < 0.25 in at least one"),
                         sliderInput("hfd_threshold", "Absolute log2FC threshold:", 
                                     min = 0, max = 3, value = 1, step = 0.1),
                         DTOutput("hfd_overlap_table"),
                         plotlyOutput("plot_hfd_4"),
                         h3("Pathways with concordant directionality between HET and MUT, and padj < 0.25 in at least one"),
                         sliderInput("hfd_gsea_threshold", "Absolute log2FC threshold:", 
                                     min = 0, max = 3, value = 1, step = 0.1),
                         DTOutput("hfd_path_overlap_table"),
                         downloadButton("download_hfd_path_table", "Download Concordant Pathways CSV"),
                         plotlyOutput("plot_hfd_5"),
                         h3("Concordance with SPONT FITC"),
                         checkboxInput("hfd_filter_concordant", "Filter for concordant direction", value = FALSE),
                         checkboxInput("hfd_apply_threshold", "Filter by effect size", value = FALSE),
                         sliderInput("hfd_threshold", "Effect size threshold", min = 0, max = 3, step = 0.1, value = 1),
                         DTOutput("hfd_comparison_table"),
                         downloadButton("hfd_download_comparison"),
                         plotlyOutput("hfd_comparison_plot")
                         
                       )
                    ),
              
              
              tabPanel("SMT_Indomethacin",
                       fluidRow(
                         column(6, plotOutput("plot_smt_indo_1")),
                         column(6, plotOutput("plot_smt_indo_2"))
                       )
              ),
              
              tabPanel("SE Supplementation",
                       fluidRow(
                         column(3, plotOutput("plot_se_supp_1")),
                         column(3, plotOutput("plot_se_supp_2")),
                         column(3, plotOutput("plot_se_supp_3")),
                         column(3, plotOutput("plot_se_supp_4"))
                       )
              ),
              
              tabPanel("SLC_ICP-MS",
                       fluidRow(
                         column(6, plotOutput("plot_slc_icp_1")),
                         column(6, plotOutput("plot_slc_icp_2"))
                       )
              )
  )
)


# ---- SERVER ----
server <- function(input, output) {
  
  full_comparisons <- list( c("WT", "HET"), c("HET","MUT"),c("WT", "MUT"))
  
  plot_slc_spont_1 <- plot_fitc(csv_filepath = here("data/phenotype/SLC_Spontaneous/FITC/SPONT_FITC_FITC_Results.csv"),
                                title_string = "FITC",
                                subtitle_string = "SPONT FITC Cohort",
                                stat_comparisons = full_comparisons) + 
    facet_wrap(~Size)
  plot_slc_spont_2 <- plot_histology(csv_filepath= here("data/phenotype/SLC_Spontaneous/FITC/SPONT_FITC_Jejunum_Histo.csv"),
                                     title_string = "Jejunum Histology",
                                     subtitle_string =  "SPONT FITC Cohort",
                                     stat_comparisons =  full_comparisons) 
  plot_slc_spont_3 <- plot_histology(csv_filepath= here("data/phenotype/SLC_ICP-MS/Ileum Spontaneous ICP-MS.csv"),
                                     title_string = "Ileum Histology",
                                     subtitle_string =  "ICP-MS Cohort",
                                     stat_comparisons =  list( c("WT", "MUT"))) 
  
  
  smt_indo <- read.csv(here("data/phenotype/SMT_Indomethacin/SMT_Indo_DAI.csv"))
  meta <- read.csv(here("data/phenotype/SMT_Indomethacin/SMT_Indo_Mouse_Info.csv")) %>%
    dplyr::select("MouseID", "Genotype","Sex")
  smt_indo <- smt_indo %>% 
    left_join(meta, by="MouseID")
  
  plot_smt_indo_1 <-   plot_avg_trajectory(smt_indo, title_string = "SMT Indomethacin DAI", subtitle_string = "LMEM Genotype p = 0.0206")
  
  
  indo <- read.csv(here("data/phenotype/SLC_Indomethacin/SLC_Indo_DAI.csv"))
  meta <- read.csv(here("data/phenotype/SLC_Indomethacin/SLC_Indomethacin_Mouse_Info.csv")) %>%
    dplyr::select("MouseID", "Genotype","Sex")
  indo <- indo %>% 
    left_join(meta, by="MouseID")
  plot_slc_indo_1 <- plot_avg_trajectory(indo, title_string = "SLC Indomethacin DAI", subtitle_string = "HET p = 0.1913, MUT p = 0.0036")
  
  plot_slc_indo_2 <- plot_histology(csv_filepath= here("data/phenotype/SLC_Indomethacin/SLC_Indo_Histo.csv"),
                                    title_string = "Ileum Histology",
                                    subtitle_string =  "SLC Indomethacin Cohort",
                                    stat_comparisons =  full_comparisons) 
  
  se_indo <- read.csv(here("data/phenotype/SE_Supp/SE_Supp_Indomethacin_DAI.csv"))
  se_indo_ctrl <- se_indo %>% filter(Diet =="Control")
  plot_se_supp_1 <- plot_avg_trajectory(se_indo_ctrl, title_string = "SE SUPP Control Diet", subtitle_string = "HET p = 0.5727, MUT p = 0.5037")
  
  se_indo_supp <- se_indo %>% filter(Diet =="Selenium")
  plot_se_supp_2 <- plot_avg_trajectory(se_indo_supp, title_string = "SE SUPP Selenium Diet", subtitle_string = "HET p = 0.7391, MUT p = 0.6525")
  
  se_histo_ctrl <- plot_histology(csv_filepath= here("data/phenotype/SE_Supp/Control_Histology.csv"),
                                  title_string = "Ileum Histology",
                                  subtitle_string =  "Indomethacin Negative Se and Ctrl Histology",
                                  stat_comparisons =  full_comparisons) +
    facet_wrap(~Diet)
  
  
  se_histo_indo <- plot_histology(csv_filepath= here("data/phenotype/SE_Supp/Se_Supp_Histology_Indomethacin_Analysis.csv"),
                                  title_string = "Ileum Histology",
                                  subtitle_string =  "Indomethacin Positive Se and Ctrl Histology",
                                  stat_comparisons =  full_comparisons) +
    facet_wrap(~Diet)
  
  # Define files and readers in a list
  files <- list(
    list(
      id = "preview01",
      path = here("results/RNA_seq/DESEQ2/DESEQ2_SPONT_FITC_MUT_vs_WT_results.csv"),
      reader = function(p) read.csv(p, row.names = 1)
    ),
    list(
      id = "preview02",
      path = here("results/RNA_seq/DESEQ2/DESEQ2_SPONT_FITC_HET_vs_WT_results.csv"),
      reader = function(p) read.csv(p, row.names = 1)
    ),
    list(
      id = "preview03",
      path = here("results/RNA_seq/GSEA/M2_GSEA_SPONT_MUT_vs_WT.csv"),
      reader = function(p) read.csv(p, row.names = 1)
    ),
    list(
      id = "preview05",
      path = here("results/RNA_seq/GSEA/M2_GSEA_SPONT_HET_vs_WT.csv"),
      reader = function(p) read.csv(p, row.names = 1)
    ),
    list(
      id = "preview07",
      path = here("results/RNA_seq/DESEQ2/DESEQ2_SMT_Neg_MUT_vs_WT_results.csv"),
      reader = function(p) read.csv(p, row.names = 1)
      ),
  
    list(
      id = "preview08",
      path = here("results/RNA_seq/GSEA/M2_GSEA_SMT_Neg_MUT_vs_WT.csv"),
      reader = function(p) read.csv(p, row.names = 1)
    ),
  list(
      id = "preview10",
      path = here("results/RNA_seq/DESEQ2/DESEQ2_HFD_MUT_vs_WT_results.csv"),
      reader = function(p) read.csv(p, row.names = 1)
      ),
  list(
    id = "preview16",
    path = here("results/RNA_seq/DESEQ2/DESEQ2_HFD_HET_vs_WT_results.csv"),
    reader = function(p) read.csv(p, row.names = 1)
  ),
  list(
    id = "preview11",
    path = here("results/RNA_seq/GSEA/M2_GSEA_HFD_Positive_MUT_vs_WT.csv"),
    reader = function(p) read.csv(p, row.names = 1)
  ),
  list(
    id = "preview17",
    path = here("results/RNA_seq/GSEA/M2_GSEA_HFD_Positive_HET_vs_WT.csv"),
    reader = function(p) read.csv(p, row.names = 1)
  ),
  list(
    id = "preview12",
    path = here("results/RNA_seq/DESEQ2/DESEQ2_STL_HET_vs_WT_results.csv"),
    reader = function(p) read.csv(p, row.names = 1)
    ),
  list(
    id = "preview13",
    path = here("results/RNA_seq/DESEQ2/DESEQ2_STL_MUT_vs_WT_results.csv"),
    reader = function(p) read.csv(p, row.names = 1)
    ),
  list(
    id = "preview14",
    path = here("results/RNA_seq/GSEA/M2_GSEA_STL_Positive_MUT_vs_WT.csv"),
    reader = function(p) read.csv(p, row.names = 1)
  ),
  list(
    id = "preview15",
    path = here("results/RNA_seq/GSEA/M2_GSEA_STL_Positive_HET_vs_WT.csv"),
    reader = function(p) read.csv(p, row.names = 1)
  )
  )
  
  # Loop over file definitions
  for (f in files) {
    local({
      id <- stri_sub(f$id, -2, -1)
      path <- f$path
      df <- f$reader(path)

      output[[paste0("filepath", id)]] <- renderText({
        paste("Loaded from:", normalizePath(path))
      })

      output[[f$id]] <- renderDT({
        datatable(df, options = list(scrollX = TRUE, pageLength = 5))
        
        #datatable(head(df, 100), options = list(scrollX = TRUE, pageLength = 5))
      })
    })
  }

  slc_spont_overlap_result <- reactive({
    make_overlap_plot(
      het_input = "results/RNA_seq/DESEQ2/DESEQ2_SPONT_FITC_HET_vs_WT_results.csv",
      mut_input = "results/RNA_seq/DESEQ2/DESEQ2_SPONT_FITC_MUT_vs_WT_results.csv",
      threshold = input$threshold,
      type = "DEG"
    )
  })
  
  slc_spont_path_overlap_result <- reactive({
    make_overlap_plot(het_input = here("results/RNA_seq/GSEA/M2_GSEA_SPONT_HET_vs_WT.csv"),
                      mut_input = here("results/RNA_seq/GSEA/M2_GSEA_SPONT_MUT_vs_WT.csv"),
                      threshold = input$gsea_threshold,
                      type="GSEA",
                      padj_cutoff=0.25)
  })
  

  
  # Render them
  output$plot_slc_spont_1 <- renderPlot({ print(plot_slc_spont_1) })
  output$plot_slc_spont_2 <- renderPlot({ print(plot_slc_spont_2 ) })
  output$plot_slc_spont_3 <- renderPlot({ print(plot_slc_spont_3 ) })

  full_wide <- reactive({
    slc_spont_overlap_result()$table
  })
  
  output$download_full_wide <- downloadHandler(
    filename = function() {
      paste0("full_wide_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(full_wide(), file, row.names = FALSE)
    }
  )
  
  path_wide <- reactive({
    slc_spont_path_overlap_result()$table
  })
  
  output$download_path_wide <- downloadHandler(
    filename = function() {
      paste0("path_wide_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(path_wide(), file, row.names = FALSE)
    }
  )
  
  output$plot_slc_spont_4 <- renderPlotly({
    slc_spont_overlap_result()$plot
  })
  output$overlap_table <- renderDT({
    datatable(slc_spont_overlap_result()$table, options = list(scrollX = TRUE, pageLength = 10))
  })
  
  output$plot_slc_spont_5<- renderPlotly({
    slc_spont_path_overlap_result()$plot
  })
  
  output$path_overlap_table <- renderDT({
    datatable(slc_spont_path_overlap_result()$table, options = list(scrollX = TRUE, pageLength = 10))
  })
  

  
  

  
  output$plot_smt_indo_1 <- renderPlot({ print(plot_smt_indo_1) })
  output$plot_slc_indo_1 <- renderPlot({ print(plot_slc_indo_1) })
  output$plot_slc_indo_2 <- renderPlot({ print(plot_slc_indo_2) })
  
  output$plot_se_supp_1 <- renderPlot({ print(plot_se_supp_1) })
  output$plot_se_supp_2 <- renderPlot({ print(plot_se_supp_2) })
  output$plot_se_supp_3 <- renderPlot({ print(se_histo_ctrl) })
  output$plot_se_supp_4 <- renderPlot({ print(se_histo_indo) })
  
 
  
  # Use the helper inside your outputs
  output$gsea_manhattan_mut <- renderPlotly({
    make_gsea_plot("results/RNA_seq/GSEA/M2_GSEA_HFD_Positive_MUT_vs_WT.csv")
  })

  output$gsea_manhattan_het <- renderPlotly({
    make_gsea_plot("results/RNA_seq/GSEA/M2_GSEA_HFD_Positive_HET_vs_WT.csv")
  })
  
  observeEvent(input$tabs, {
    if (input$tabs == "SMT_Negative") {
      plot_smt_negative_1 <- plot_histology(csv_filepath= here("data/phenotype/SMT_Indomethacin/SMT_Both_Histo.csv"),
                                            title_string = "Ileum Histology",
                                            subtitle_string =  "SMT Cohort",
                                            stat_comparisons =  list( c("WT", "MUT"))) +
        facet_wrap(~Treatment)
      plot_smt_negative_2 <- plot_fitc(csv_filepath = here("data/phenotype/SMT_Negative/SMT_Neg_FITC.csv"),
                                       title_string = "FITC",
                                       subtitle_string = "SMT Negative FITC",
                                       stat_comparisons =  list( c("WT", "MUT"))) + 
        facet_wrap(~Size)
      
      # plot_smt_negative_3 <- readRDS(here("results/SMT_Shotgun_DAT.RDS"))
      # plot_smt_negative_4 <- readRDS(here("results/SMT_Shotgun_Beta_Diversity.RDS"))
      # plot_smt_negative_5 <- readRDS(here("results/phenotype/SMT_Cell_Frac.RDS"))
      #table_smt_negative <- read.csv(here("results/phenotype/SMT_Neg_Cell_Frac_Wilcox.csv"), row.names=1)
      
      output$plot_smt_negative_1 <- renderPlot({ print(plot_smt_negative_1 ) })
      output$plot_smt_negative_2 <- renderPlot({ print(plot_smt_negative_2 ) })
      # output$plot_smt_negative_3 <- renderPlot({ print(plot_smt_negative_3 ) })
      # output$plot_smt_negative_4 <- renderPlot({ print(plot_smt_negative_4 ) })
      # output$plot_smt_negative_5 <- renderPlot({ print(plot_smt_negative_5 ) })
      #output$SMT_Neg_wilcox_table <- renderTable(table_smt_negative)
      
      output$SMT_Neg_gsea_manhattan <- renderPlotly({
        make_gsea_plot("results/RNA_seq/GSEA/M2_GSEA_SMT_Neg_MUT_vs_WT.csv")
      })
      
      comparison_results <- reactive({
        make_comparison_plot_smt(
          gsea_file_smt = here("results/RNA_seq/GSEA/M2_GSEA_SMT_Neg_MUT_vs_WT.csv"),
          gsea_file_spont = here("results/RNA_seq/GSEA/SPONT_FITC_PATH_CONCORDANT.csv"),
          threshold = input$threshold,
          filter_concordant = input$filter_concordant,
          apply_threshold = input$apply_threshold
        )
      })
      
      output$comparison_plot <- renderPlotly({
        comparison_results()$plot
      })
      
      output$comparison_table <- renderDT({
        datatable(comparison_results()$table, options = list(scrollX = TRUE))
      })
      
      
      output$download_comparison <- downloadHandler(
        filename = function() paste0("comparison_", Sys.Date(), ".csv"),
        content = function(file) write.csv(comparison_result()$full_wide, file, row.names = FALSE)
      )

    }
  })
  
  observeEvent(input$tabs, {
    if (input$tabs == "SLC_TL1A") {
      plot_tl1a_1 <- plot_histology(csv_filepath= here("data/phenotype/SLC_TL1A/STL_Histology_Ileum.csv"),
                                    title_string = "Ileum Histology",
                                    subtitle_string =  "STL Cohort",
                                    stat_comparisons =  full_comparisons) +
        facet_wrap(~TL1A)
      
      tl1a_overlap_result <- reactive({
        make_overlap_plot(
          het_input = "results/RNA_seq/DESEQ2/DESEQ2_STL_HET_vs_WT_results.csv",
          mut_input = "results/RNA_seq/DESEQ2/DESEQ2_STL_MUT_vs_WT_results.csv",
          threshold = input$threshold,
          type = "DEG"
        )
      })
      
      tl1a_path_overlap_result <- reactive({
        make_overlap_plot(het_input = here("results/RNA_seq/GSEA/M2_GSEA_STL_Positive_HET_vs_WT.csv"),
                          mut_input = here("results/RNA_seq/GSEA/M2_GSEA_STL_Positive_MUT_vs_WT.csv"),
                          threshold = input$gsea_threshold,
                          type="GSEA",
                          padj_cutoff=0.25)
      })
      
  
      output$plot_tl1a_1 <- renderPlot({ print(plot_tl1a_1) })
      output$plot_tl1a_4 <- renderPlotly({
        tl1a_overlap_result()$plot
      })
      output$tl1a_overlap_table <- renderDT({
        datatable(tl1a_overlap_result()$table, options = list(scrollX = TRUE, pageLength = 10))
      })
      
      output$plot_tl1a_5<- renderPlotly({
        tl1a_path_overlap_result()$plot
      })
      
      output$tl1a_path_overlap_table <- renderDT({
        datatable(tl1a_path_overlap_result()$table, options = list(scrollX = TRUE, pageLength = 10))
      })
      
      output$download_tl1a_path_table <- downloadHandler(
        filename = function() {
          paste0("tl1a_path_table_", Sys.Date(), ".csv")
        },
        content = function(file) {
          write.csv(tl1a_path_overlap_result()$table, file, row.names = FALSE)
        }
      )
      
      tl1a_comparison_results <- reactive({
        make_comparison_plot(
          #gsea_file_1 = here("results/RNA_seq/GSEA/HFD_PATH_CONCORDANT.csv"),
          gsea_file_1 = tl1a_path_overlap_result()$table,
          gsea_file_spont = here("results/RNA_seq/GSEA/SPONT_FITC_PATH_CONCORDANT.csv"),
          threshold = input$tl1a_threshold,
          filter_concordant = input$tl1a_filter_concordant,
          apply_threshold = input$tl1a_apply_threshold
        )
      })
      
      output$tl1a_comparison_plot <- renderPlotly({
        tl1a_comparison_results()$plot
      })
      
      output$tl1a_comparison_table <- renderDT({
        datatable(tl1a_comparison_results()$table, options = list(scrollX = TRUE))
      })
      
      
      output$tl1a_download_comparison <- downloadHandler(
        filename = function() paste0("comparison_", Sys.Date(), ".csv"),
        content = function(file) write.csv(tl1a_comparison_results()$full_wide, file, row.names = FALSE)
      )
      
    }
  })
  
  observeEvent(input$tabs, {
    if (input$tabs == "SLC_HFD") {
      output$plot_slc_hfd_1 <- renderPlot({ print(plot_hfd_1) })
      output$plot_slc_hfd_2 <- renderPlot({ print(plot_hfd_2) })
      plot_hfd_1 <- plot_histology(csv_filepath= here("data/phenotype/SLC_HFD/HFD_Histo.csv"),
                                   title_string = "Ileum Histology",
                                   subtitle_string =  "HFD Cohort",
                                   stat_comparisons =  full_comparisons) 
      plot_hfd_2 <- plot_fitc(csv_filepath = here("data/phenotype/SLC_HFD/HFD_FITC_Analysis.csv"),
                              title_string = "FITC",
                              subtitle_string = "HFD Cohort",
                              stat_comparisons = full_comparisons) + 
        facet_wrap(~Size)
      
      hfd_overlap_result <- reactive({
        make_overlap_plot(
          het_input = "results/RNA_seq/DESEQ2/DESEQ2_HFD_HET_vs_WT_results.csv",
          mut_input = "results/RNA_seq/DESEQ2/DESEQ2_HFD_MUT_vs_WT_results.csv",
          threshold = input$hfd_threshold,
          type = "DEG"
        )
      })
      
      hfd_path_overlap_result <- reactive({
        make_overlap_plot(het_input = here("results/RNA_seq/GSEA/M2_GSEA_HFD_Positive_HET_vs_WT.csv"),
                          mut_input = here("results/RNA_seq/GSEA/M2_GSEA_HFD_Positive_MUT_vs_WT.csv"),
                          threshold = input$hfd_gsea_threshold,
                          type="GSEA",
                          padj_cutoff=0.25)
      })
      
      output$plot_hfd_4 <- renderPlotly({
        hfd_overlap_result()$plot
      })
      output$hfd_overlap_table <- renderDT({
        datatable(hfd_overlap_result()$table, options = list(scrollX = TRUE, pageLength = 10))
      })
      
      output$plot_hfd_5<- renderPlotly({
        hfd_path_overlap_result()$plot
      })
      
      output$hfd_path_overlap_table <- renderDT({
        datatable(hfd_path_overlap_result()$table, options = list(scrollX = TRUE, pageLength = 10))
      })
      
      output$download_hfd_path_table <- downloadHandler(
        filename = function() {
          paste0("hfd_path_table_", Sys.Date(), ".csv")
        },
        content = function(file) {
          write.csv(hfd_path_overlap_result()$table, file, row.names = FALSE)
        }
      )
      
      hfd_comparison_results <- reactive({
        make_comparison_plot(
          #gsea_file_1 = here("results/RNA_seq/GSEA/HFD_PATH_CONCORDANT.csv"),
          gsea_file_1 = hfd_path_overlap_result()$table,
          gsea_file_spont = here("results/RNA_seq/GSEA/SPONT_FITC_PATH_CONCORDANT.csv"),
          threshold = input$hfd_threshold,
          filter_concordant = input$hfd_filter_concordant,
          apply_threshold = input$hfd_apply_threshold
        )
      })

      output$hfd_comparison_plot <- renderPlotly({
        hfd_comparison_results()$plot
      })

      output$hfd_comparison_table <- renderDT({
        datatable(hfd_comparison_results()$table, options = list(scrollX = TRUE))
      })


      output$hfd_download_comparison <- downloadHandler(
        filename = function() paste0("comparison_", Sys.Date(), ".csv"),
        content = function(file) write.csv(hfd_comparison_results()$full_wide, file, row.names = FALSE)
      )
      
    }
  })
  
  
}

# ---- RUN APP ----
shinyApp(ui, server)
