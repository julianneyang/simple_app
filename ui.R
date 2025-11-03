

source("global.R")

# ---- UI ----
ui <- fluidPage(
  titlePanel("Aggregated Data from the SLC SI Project"),
  tabsetPanel(id = "tabs",
              tabPanel("SLC_Spontaneous",
                       fluidRow(
                         h2("Phenotype Results"),
                         column(3, plotOutput("plot_slc_spont_1A")),
                         column(3, plotOutput("plot_slc_spont_1B")),
                         column(3, plotOutput("plot_slc_spont_2")),
                         column(3, plotOutput("plot_slc_spont_3"))
                       ),
                       fluidRow(
                         column(3, h3("150 kDa Sample Size", align = "center"), tableOutput("spont_fitc_summary_table")),
                         column(3, h3("4 KDa Sample Size"), tableOutput("spont_fitc_4_summary_table")),
                         column(3, plotOutput("plot_slc_spont_1C")),
                         column(3, plotOutput("plot_slc_spont_1D"))
                       ),
                       fluidRow(
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
                         plotlyOutput("plot_slc_spont_5"),
                         plotOutput("plot_slc_spont_proteomics")
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
                         column(3, plotOutput("plot_se_supp_4")),
                       ),
                       fluidRow(
                         column(4, plotOutput("plot_se_supp_7")),
                         column(4, plotOutput("plot_se_supp_8")),
                         column(4, plotOutput("plot_se_supp_10")),
                       ),
                       h4("Combining Control Diet with SPONT FITC"),
                       fluidRow(
                         column(4, plotOutput("plot_se_supp_11")),
                         column(8, plotOutput("plot_se_supp_12"))
                       ),
                       h4("Wilcoxon rank sum comparisons due to outliers"),
                       fluidRow(
                         
                         column(12, plotOutput("plot_se_supp_5")),
                       ),
                       fluidRow(
                         column(12, plotOutput("plot_se_supp_6")),
                       ),
                       fluidRow(
                         column(12, plotOutput("plot_se_supp_9")),
                       )
              ),
              
              tabPanel("SLC_ICP-MS",
                       h4("Wilcoxon rank sum comparisons due to outliers"),
                       fluidRow(
                         column(12, plotOutput("plot_icp_1"))
                        )
              ),
              tabPanel(
                "SPONT_Proteomics",         
                selectInput(
                  inputId = "proteomics_folder",
                  label = "Choose dataset:",
                  choices = c("2025_only", "Combined"),
                  selected = "2025_only"
                ),
                uiOutput("tables_ui")
              )
  )
)





