renv::restore()
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

# ---- Load data ----
here::i_am("app.R")
# ---- UI ----
ui <- fluidPage(
  titlePanel("Aggregated Data from the SLC SI Project"),
  
  tabsetPanel(
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
               h3("GSEA: MUT vs WT, M7 Immunologic Signature Gene Sets"),
               textOutput("filepath4"),
               DTOutput("preview4"),
               plotlyOutput("gsea_manhattan", height = "600px"),
               h3("GSEA: HET vs WT"),
               textOutput("filepath5"),
               DTOutput("preview5"),
               h3("GSEA: HET vs WT, M7 Immunologic Signature Gene Sets"),
               textOutput("filepath6"),
               DTOutput("preview6"),
               plotlyOutput("gsea_manhattan_het", height = "600px"),
               
             )
             
    ),
    tabPanel("SMT_Negative",
             fluidRow(
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
               h3("GSEA: MUT vs WT, M7 Immunologic Signature Gene Sets"),
               textOutput("filepath9"),
               DTOutput("preview9"),
               plotlyOutput("SMT_Neg_gsea_manhattan", height = "600px"),
             )
    ),
    
    tabPanel("SLC_TL1A",
             fluidRow(
               column(6, plotOutput("plot_slc_tl1a_1")),
               column(6, plotOutput("plot_slc_tl1a_2"))
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
               column(6, plotOutput("plot_slc_hfd_1")),
               column(6, plotOutput("plot_slc_hfd_2"))
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
  
  # Read in your pre-made ggplot objects
  plot_slc_spont_1 <- readRDS(here("results/phenotype/Spont_FITC_4KDa.RDS"))
  plot_slc_spont_2 <- readRDS(here("results/phenotype/Spont_FITC_Jej_Histo.RDS"))
  plot_slc_spont_3 <- readRDS(here("results/phenotype/ICP_MS_Histo.RDS"))
  plot_smt_negative_1 <- readRDS(here("results/phenotype/SMT_Ileum_Histo.RDS"))
  plot_smt_negative_2 <- readRDS(here("results/phenotype/SMT_Neg_FITC.RDS"))
  plot_smt_negative_3 <- readRDS(here("results/SMT_Shotgun_DAT.RDS"))
  plot_smt_negative_4 <- readRDS(here("results/SMT_Shotgun_Beta_Diversity.RDS"))
  
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
    list(
      id = "preview3",
      path = here("results/RNA_seq/GSEA/GSEA_SPONT_FITC_MUT_vs_WT.RDS"),
      reader = function(p) readRDS(p) %>% arrange(padj)
    ),
    list(
      id = "preview4",
      path = here("results/RNA_seq/GSEA/M7_GSEA_SPONT_FITC_MUT_vs_WT.RDS"),
      reader = function(p) readRDS(p) %>% arrange(padj)
    ),
    list(
      id = "preview5",
      path = here("results/RNA_seq/GSEA/GSEA_SPONT_FITC_HET_vs_WT.RDS"),
      reader = function(p) readRDS(p) %>% arrange(padj)
    ),
    list(
      id = "preview6",
      path = here("results/RNA_seq/GSEA/M7_GSEA_SPONT_FITC_HET_vs_WT.RDS"),
      reader = function(p) readRDS(p) %>% arrange(padj)
    ),
    list(
      id = "preview7",
      path = here("results/RNA_seq/DESEQ2/DESEQ2_SMT_Neg_MUT_vs_WT_results.csv"),
      reader = function(p) read.csv(p, row.names = 1)
    ),
    list(
      id = "preview8",
      path = here("results/RNA_seq/GSEA/GSEA_SMT_Neg_MUT_vs_WT.RDS"),
      reader = function(p) readRDS(p) %>% arrange(padj)
    ),
    list(
      id = "preview9",
      path = here("results/RNA_seq/GSEA/M7_GSEA_SMT_Neg_MUT_vs_WT.RDS"),
      reader = function(p) readRDS(p) %>% arrange(padj)
    )
  )
  
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
        datatable(df, options = list(scrollX = TRUE, pageLength = 5))
      })
    })
  }
  
  
  # Render them
  output$plot_slc_spont_1 <- renderPlot({ print(plot_slc_spont_1) })
  output$plot_slc_spont_2 <- renderPlot({ print(plot_slc_spont_2 ) })
  output$plot_slc_spont_3 <- renderPlot({ print(plot_slc_spont_3 ) })
  output$plot_smt_negative_1 <- renderPlot({ print(plot_smt_negative_1 ) })
  output$plot_smt_negative_2 <- renderPlot({ print(plot_smt_negative_2 ) })
  output$plot_smt_negative_3 <- renderPlot({ print(plot_smt_negative_3 ) })
  output$plot_smt_negative_4 <- renderPlot({ print(plot_smt_negative_4 ) })
  
  output$plot_tl1a2 <- renderPlot({ print(plot_tl1a2) })
  
  set.seed(123)
  gsea_df <- data.frame(readRDS(here("results/RNA_seq/GSEA/M7_GSEA_SPONT_FITC_MUT_vs_WT.RDS"))) %>%
    mutate(logFDR = -log10(padj))
  
  output$gsea_manhattan <- renderPlotly({
    plot_ly(
      data = gsea_df,
      x = ~NES,
      y = ~logFDR,
      text = ~paste0(
        "Pathway: ", pathway, "<br>",
        "NES: ", round(NES, 2), "<br>",
        "FDR: ", signif(padj, 3)
      ),
      type = "scatter",
      mode = "markers",
      marker = list(size = 10,
                    color = ~logFDR,
                    colorscale = "Viridis",
                    showscale = TRUE)
    ) %>%
      layout(
        xaxis = list(title = "Normalized Enrichment Score (NES)"),
        yaxis = list(title = "-log10(padj)"),
        hovermode = "closest"
      )
  })
  
  gsea_df_het <- data.frame(readRDS(here("results/RNA_seq/GSEA/M7_GSEA_SPONT_FITC_HET_vs_WT.RDS"))) %>%
    mutate(logFDR = -log10(padj))
  
  output$gsea_manhattan_het <- renderPlotly({
    plot_ly(
      data = gsea_df_het,
      x = ~NES,
      y = ~logFDR,
      text = ~paste0(
        "Pathway: ", pathway, "<br>",
        "NES: ", round(NES, 2), "<br>",
        "FDR: ", signif(padj, 3)
      ),
      type = "scatter",
      mode = "markers",
      marker = list(size = 10,
                    color = ~logFDR,
                    colorscale = "Viridis",
                    showscale = TRUE)
    ) %>%
      layout(
        xaxis = list(title = "Normalized Enrichment Score (NES)"),
        yaxis = list(title = "-log10(padj)"),
        hovermode = "closest"
      )
  })
  
  gsea_df_smt <- data.frame(readRDS(here("results/RNA_seq/GSEA/M7_GSEA_SMT_Neg_MUT_vs_WT.RDS"))) %>%
    mutate(logFDR = -log10(padj))
  
  output$SMT_Neg_gsea_manhattan <- renderPlotly({
    plot_ly(
      data = gsea_df_smt,
      x = ~NES,
      y = ~logFDR,
      text = ~paste0(
        "Pathway: ", pathway, "<br>",
        "NES: ", round(NES, 2), "<br>",
        "FDR: ", signif(padj, 3)
      ),
      type = "scatter",
      mode = "markers",
      marker = list(size = 10,
                    color = ~logFDR,
                    colorscale = "Viridis",
                    showscale = TRUE)
    ) %>%
      layout(
        xaxis = list(title = "Normalized Enrichment Score (NES)"),
        yaxis = list(title = "-log10(padj)"),
        hovermode = "closest"
      )
  })
  
}

# ---- RUN APP ----
shinyApp(ui, server)
