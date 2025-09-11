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
                         textOutput("filepath1"),
                         DTOutput("preview1"),
                         h3("DESeq2: HET vs WT"),
                         textOutput("filepath2"),
                         DTOutput("preview2"),
                         br(),
                         h3("GSEA: MUT vs WT"),
                         textOutput("filepath3"),
                         DTOutput("preview3"),
                         plotlyOutput("gsea_manhattan_mut", height = "600px"),
                         h3("GSEA: MUT vs WT, M7 Immunologic Signature Gene Sets"),
                         textOutput("filepath4"),
                         DTOutput("preview4"),
                         plotlyOutput("M7_gsea_manhattan_mut", height = "600px"),
                         h3("GSEA: HET vs WT"),
                         textOutput("filepath5"),
                         DTOutput("preview5"),
                         plotlyOutput("gsea_manhattan_het", height = "600px"),
                         h3("GSEA: HET vs WT, M7 Immunologic Signature Gene Sets"),
                         textOutput("filepath6"),
                         DTOutput("preview6"),
                         plotlyOutput("M7_gsea_manhattan_het", height = "600px"),
                         br(),
                         plotOutput("plot_slc_spont_4")
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
                         textOutput("filepath7"),
                         DTOutput("preview7"),
                         br(),
                         h3("GSEA: MUT vs WT"),
                         textOutput("filepath8"),
                         DTOutput("preview8"),
                         plotlyOutput("SMT_Neg_gsea_manhattan", height = "600px"),
                         h3("GSEA: MUT vs WT, M7 Immunologic Signature Gene Sets"),
                         textOutput("filepath9"),
                         DTOutput("preview9"),
                         plotlyOutput("M7_SMT_Neg_gsea_manhattan", height = "600px"),
                         br(),
                         tableOutput("SMT_Neg_wilcox_table"),
                         plotOutput("plot_smt_negative_5")
                       )
              ),
              
              tabPanel("SLC_TL1A",
                       fluidRow(
                         column(6, plotOutput("plot_tl1a_1")),
                         column(6, plotOutput("plot_tl1a_2")),
                         h2("RNA-sequencing Results"),
                         br(),
                         h3("DESeq2: HET vs WT"),
                         textOutput("filepath11"),
                         DTOutput("preview11"),
                         h3("DESeq2: MUT vs WT"),
                         textOutput("filepath12"),
                         DTOutput("preview12")
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
                       )
              ),
              
              tabPanel("SMT_Indomethacin",
                       fluidRow(
                         column(6, plotOutput("plot_smt_indo_1")),
                         column(6, plotOutput("plot_smt_indo_2"))
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
  
  # Read in your pre-made ggplot objects
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
  
  #plot_spont_cibersort <- readRDS(here("results/phenotype/SPONT_Cell_Frac.RDS"))
  plot_tl1a_1 <- plot_histology(csv_filepath= here("data/phenotype/SLC_TL1A/STL_Histology_Ileum.csv"),
                                title_string = "Ileum Histology",
                                subtitle_string =  "STL Cohort",
                                stat_comparisons =  full_comparisons) +
    facet_wrap(~TL1A)
  
  plot_hfd_1 <- plot_histology(csv_filepath= here("data/phenotype/SLC_HFD/HFD_Histo.csv"),
                               title_string = "Ileum Histology",
                               subtitle_string =  "HFD Cohort",
                               stat_comparisons =  full_comparisons) 
  plot_hfd_2 <- plot_fitc(csv_filepath = here("data/phenotype/SLC_HFD/HFD_FITC_Analysis.csv"),
                          title_string = "FITC",
                          subtitle_string = "HFD Cohort",
                          stat_comparisons = full_comparisons) + 
    facet_wrap(~Size)
  
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
  
  
  # Define files and readers in a list
  files <- list(
    list(
      id = "preview1",
      path = here("results/RNA_seq/DESEQ2/DESEQ2_SPONT_FITC_MUT_vs_WT_results.csv"),
      reader = function(p) read.csv(p, row.names = 1)
    ),
    list(
      id = "preview2",
      path = here("results/RNA_seq/DESEQ2/DESEQ2_SPONT_FITC_HET_vs_WT_results.csv"),
      reader = function(p) read.csv(p, row.names = 1)
    ),
    # list(
    #   id = "preview3",
    #   path = here("results/RNA_seq/GSEA/GSEA_SPONT_FITC_MUT_vs_WT.RDS"),
    #   reader = function(p) readRDS(p) %>% arrange(padj)
    # ),
    # list(
    #   id = "preview4",
    #   path = here("results/RNA_seq/GSEA/M7_GSEA_SPONT_FITC_MUT_vs_WT.RDS"),
    #   reader = function(p) readRDS(p) %>% arrange(padj)
    # ),
    # list(
    #   id = "preview5",
    #   path = here("results/RNA_seq/GSEA/GSEA_SPONT_FITC_HET_vs_WT.RDS"),
    #   reader = function(p) readRDS(p) %>% arrange(padj)
    # ),
    # list(
    #   id = "preview6",
    #   path = here("results/RNA_seq/GSEA/M7_GSEA_SPONT_FITC_HET_vs_WT.RDS"),
    #   reader = function(p) readRDS(p) %>% arrange(padj)
    # ),
    list(
      id = "preview7",
      path = here("results/RNA_seq/DESEQ2/DESEQ2_SMT_Neg_MUT_vs_WT_results.csv"),
      reader = function(p) read.csv(p, row.names = 1)
      ),
  
  #   list(
  #     id = "preview8",
  #     path = here("results/RNA_seq/GSEA/GSEA_SMT_Neg_MUT_vs_WT.RDS"),
  #     reader = function(p) readRDS(p) %>% arrange(padj)
  #   ),
  #   list(
  #     id = "preview9",
  #     path = here("results/RNA_seq/GSEA/M7_GSEA_SMT_Neg_MUT_vs_WT.RDS"),
  #     reader = function(p) readRDS(p) %>% arrange(padj)
  #   )
  list(
      id = "preview10",
      path = here("results/RNA_seq/DESEQ2/DESEQ2_HFD_MUT_vs_WT_results.csv"),
      reader = function(p) read.csv(p, row.names = 1)
      ),
  list(
    id = "preview11",
    path = here("results/RNA_seq/DESEQ2/DESEQ2_STL_Positive_HET_vs_WT_results.csv"),
    reader = function(p) read.csv(p, row.names = 1)
    ),
  list(
    id = "preview12",
    path = here("results/RNA_seq/DESEQ2/DESEQ2_STL_Positive_MUT_vs_WT_results.csv"),
    reader = function(p) read.csv(p, row.names = 1)
    )
  )
  # 
  # Loop over file definitions
  for (f in files) {
    local({
      id <- stri_sub(f$id, -1, -1)
      path <- f$path
      df <- f$reader(path)

      output[[paste0("filepath", id)]] <- renderText({
        paste("Loaded from:", normalizePath(path))
      })

      output[[f$id]] <- renderDT({
        datatable(head(df, 100), options = list(scrollX = TRUE, pageLength = 5))
      })
    })
  }

  
  # Render them
  output$plot_slc_spont_1 <- renderPlot({ print(plot_slc_spont_1) })
  output$plot_slc_spont_2 <- renderPlot({ print(plot_slc_spont_2 ) })
  output$plot_slc_spont_3 <- renderPlot({ print(plot_slc_spont_3 ) })
  #output$plot_slc_spont_4 <- renderPlot({ print(plot_spont_cibersort ) })
  
  output$plot_tl1a_1 <- renderPlot({ print(plot_tl1a_1) })
  output$plot_slc_hfd_1 <- renderPlot({ print(plot_hfd_1) })
  output$plot_slc_hfd_2 <- renderPlot({ print(plot_hfd_2) })
  
  output$plot_smt_indo_1 <- renderPlot({ print(plot_smt_indo_1) })
  output$plot_slc_indo_1 <- renderPlot({ print(plot_slc_indo_1) })
  
  # set.seed(123)
  # 
  # # Helper function: read, preprocess, and make GSEA Manhattan plot
  # make_gsea_plot <- function(rds_path) {
  #   df <- readRDS(here(rds_path)) %>%
  #     as.data.frame() %>%
  #     mutate(logFDR = -log10(padj))
  #   
  #   plot_ly(
  #     data = df,
  #     x = ~NES,
  #     y = ~logFDR,
  #     text = ~paste0(
  #       "Pathway: ", pathway, "<br>",
  #       "NES: ", round(NES, 2), "<br>",
  #       "FDR: ", signif(padj, 3)
  #     ),
  #     type = "scatter",
  #     mode = "markers",
  #     marker = list(
  #       size = 10,
  #       color = ~logFDR,
  #       colorscale = "Viridis",
  #       showscale = TRUE
  #     )
  #   ) %>%
  #     layout(
  #       xaxis = list(title = "Normalized Enrichment Score (NES)"),
  #       yaxis = list(title = "-log10(padj)"),
  #       hovermode = "closest"
  #     )
  # }
  # 
  # # Use the helper inside your outputs
  # output$M7_gsea_manhattan_mut <- renderPlotly({
  #   make_gsea_plot("results/RNA_seq/GSEA/M7_GSEA_SPONT_FITC_MUT_vs_WT.RDS")
  # })
  # 
  # output$M7_gsea_manhattan_het <- renderPlotly({
  #   make_gsea_plot("results/RNA_seq/GSEA/M7_GSEA_SPONT_FITC_HET_vs_WT.RDS")
  # })
  # 
  # output$gsea_manhattan_mut <- renderPlotly({
  #   make_gsea_plot("results/RNA_seq/GSEA/GSEA_SPONT_FITC_MUT_vs_WT.RDS")
  # })
  # 
  # output$gsea_manhattan_het <- renderPlotly({
  #   make_gsea_plot("results/RNA_seq/GSEA/GSEA_SPONT_FITC_HET_vs_WT.RDS")
  # })
  
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
      output$SMT_Neg_wilcox_table <- renderTable(table_smt_negative)
      
      # output$SMT_Neg_gsea_manhattan <- renderPlotly({
      #   make_gsea_plot("results/RNA_seq/GSEA/GSEA_SMT_Neg_MUT_vs_WT.RDS")
      # })
      # 
      # output$M7_SMT_Neg_gsea_manhattan <- renderPlotly({
      #   make_gsea_plot("results/RNA_seq/GSEA/M7_GSEA_SMT_Neg_MUT_vs_WT.RDS")
      # })
    }
  })
  
}

# ---- RUN APP ----
shinyApp(ui, server)
