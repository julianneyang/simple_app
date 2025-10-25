
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
library(pheatmap)

here::i_am("global.R")

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

plot_fitc <- function(dataframe, title_string, subtitle_string, stat_comparisons = list( c("WT", "HET"), c("HET","MUT"),c("WT", "MUT")),
                      color_vector = c("WT"="black","HET"= "navy","MUT"="firebrick")) {
  df <- as.data.frame(dataframe)
  df$Genotype <- factor(df$Genotype, levels=c("WT", "HET", "MUT"))
  
  df_150 <- df %>% filter(Size=="150_kDa")
  df_4 <- df %>% filter(Size=="4_kDa")
  
  if(dim(df_150)[1]!=0){
    lm <- lm(Plasma_FITC ~ Sex + Genotype, data = df_150)
    print("150 kDa")
    print(summary(lm))
  }
  
  if(dim(df_4)[1]!=0){
    lm_4kDa <- lm(Plasma_FITC ~ Sex + Genotype, data = df_4)
    print("4 kDa")
    print(summary(lm_4kDa))
  }
  
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

plot_fitc_diet_comparison <- function(dataframe, title_string, subtitle_string, stat_comparisons = list( c("WT", "HET"), c("HET","MUT"),c("WT", "MUT")),
                      color_vector = c("WT"="black","HET"= "navy","MUT"="firebrick")) {
  df <- as.data.frame(dataframe)
  df$Genotype <- factor(df$Genotype, levels=c("WT", "HET", "MUT"))
  
  df_WT <- df %>% filter(Genotype=="WT")
  df_HET <- df %>% filter(Size=="HET")
  df_MUT <- df %>% filter(Size=="MUT")
  
  if(dim(df_WT)[1]!=0){
    lm <- lm(Plasma_FITC ~ Sex + Diet, data = df_WT)
    print("WT Se vs Ctrl")
    print(summary(lm))
  }
  
  if(dim(df_HET)[1]!=0){
    lm <- lm(Plasma_FITC ~ Sex + Diet, data = df_HET)
    print("HET Se vs Ctrl")
    print(summary(lm))
  }
  
  if(dim(df_MUT)[1]!=0){
    lm <- lm(Plasma_FITC ~ Sex + Diet, data = df_MUT)
    print("MUT Se vs Ctrl")
    print(summary(lm))
  }
  
  
  ggplot(df, aes(x = Diet_Genotype, y = Plasma_FITC, fill = Genotype)) +
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



set.seed(123)
make_icp_ms_plot <- function(df, subtitle_string, 
                             stat_comparisons = list( c("WT", "HET"), c("HET","MUT"),c("WT", "MUT")),
                             color_vector = c("WT"="black","HET"= "navy","MUT"="firebrick")){
  
  df <- as.data.frame(df)
  df$Genotype <- factor(df$Genotype, levels=c("WT","HET","MUT"))
  
  # Split into list by element
  element_list <- split(df, df$Element)
  
  
  
  plots <- lapply(seq_along(element_list), function(i) {
    maximum <- max(element_list[[i]]$Concentration)
    print(summary(element_list[[i]]$Concentration))
    new_ylim <- c(0, 1.5*maximum)
    print(new_ylim)
    
    ggplot(element_list[[i]], aes(x = Genotype, y = Concentration, fill = Genotype)) +
      geom_boxplot(alpha=0.5) +
      scale_fill_manual(values=color_vector) +
      geom_jitter(width = 0.05, alpha=0.8) +
      theme_cowplot(12) +
      theme(legend.position = "none")+
      ggtitle(as.character(element_list[[i]]$Element))+
      labs(subtitle = {{subtitle_string}})+
      stat_compare_means(comparisons = {{stat_comparisons}},method = "wilcox",
                         step.increase = 0.2)+
      ylab("Concentration [ug/g]")+
      xlab(NULL)+
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust=0.5)) + 
      scale_y_continuous(limits = new_ylim) 
    
  })
  
  agg_fig <- plot_grid(plotlist = plots, nrow = 1,  ncol = 7)
  return(agg_fig)
}

make_icp_ms_plot_diet_comparison <- function(df, subtitle_string, 
                                             stat_comparisons = list( c("Control_WT", "Selenium_WT"), 
                                                                      c("Control_HET","Selenium_HET"),
                                                                      c("Control_MUT", "Selenium_MUT")),
                                             color_vector = c("WT"="black","HET"= "navy","MUT"="firebrick")){
  
  df <- as.data.frame(df)
  df$Diet_Genotype <- paste0(df$Diet, "_", df$Genotype)
  df$Diet_Genotype <- factor(df$Diet_Genotype, levels=c(
    "Control_WT", "Selenium_WT",
    "Control_HET","Selenium_HET",
    "Control_MUT", "Selenium_MUT"))
  
  # Split into list by element
  element_list <- split(df, df$Element)
  
  
  
  plots <- lapply(seq_along(element_list), function(i) {
    maximum <- max(element_list[[i]]$Concentration)
    print(summary(element_list[[i]]$Concentration))
    new_ylim <- c(0, 1.5*maximum)
    print(new_ylim)
    
    ggplot(element_list[[i]], aes(x = Diet_Genotype, y = Concentration, fill = Genotype)) +
      geom_boxplot(alpha=0.5) +
      scale_fill_manual(values=color_vector) +
      geom_jitter(width = 0.05, alpha=0.8) +
      theme_cowplot(12) +
      theme(legend.position = "none") + 
      ggtitle(as.character(element_list[[i]]$Element))+
      labs(subtitle = {{subtitle_string}})+
      stat_compare_means(comparisons = {{stat_comparisons}},method = "wilcox",
                         step.increase = 0.2)+
      ylab("Concentration [ug/g]")+
      xlab(NULL)+
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust=0.5)) + 
      scale_y_continuous(limits = new_ylim) +
      theme(axis.text.x = element_text(angle = 60, vjust=0.1))
    
  })
  
  agg_fig <- plot_grid(plotlist = plots, nrow = 1,  ncol = 7)
  return(agg_fig)
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