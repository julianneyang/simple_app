

source("global.R")


# ---- SERVER ----
server <- function(input, output) {
  
  full_comparisons <- list( c("WT", "HET"), c("HET","MUT"),c("WT", "MUT"))
  
  spont_fitc_df <- read.csv(here("data/phenotype/SLC_Spontaneous/FITC/SPONT_FITC_FITC_Results.csv"))
  plot_slc_spont_1A <- plot_fitc(dataframe = spont_fitc_df %>% filter(Size =="150_kDa"),
                                title_string = "FITC 150 kDa",
                                subtitle_string = "SPONT FITC Cohort",
                                stat_comparisons = full_comparisons) 
  
  plot_slc_spont_1C <- plot_fitc(dataframe = spont_fitc_df %>% filter(Size =="150_kDa"),
                                 title_string = "FITC 150 kDa",
                                 subtitle_string = "SPONT FITC Cohort",
                                 stat_comparisons = full_comparisons) + 
    facet_wrap(~Sex)
  
  # Summarize counts by Sex and Genotype
  summary_fitc <- spont_fitc_df %>% 
    filter(Size =="150_kDa") %>%
    group_by(Sex, Genotype) %>%
    summarise(Num_Mice = n_distinct(MouseID), .groups = "drop")
  
  # Render as table
  output$spont_fitc_summary_table <- renderTable({
    summary_fitc
  })
  
  plot_slc_spont_1B <- plot_fitc(dataframe = spont_fitc_df %>% filter(Size=="4_kDa"),
                                title_string = "FITC",
                                subtitle_string = "SPONT FITC Cohort",
                                stat_comparisons = full_comparisons) + 
    facet_wrap(~Size)
  
  plot_slc_spont_1D <- plot_fitc(dataframe = spont_fitc_df %>% filter(Size=="4_kDa"),
                                 title_string = "FITC",
                                 subtitle_string = "SPONT FITC Cohort",
                                 stat_comparisons = full_comparisons) + 
    facet_wrap(~Sex + Size)
  
  # Summarize counts by Sex and Genotype
  summary_fitc_4 <- spont_fitc_df %>% 
    filter(Size =="4_kDa") %>%
    group_by(Sex, Genotype) %>%
    summarise(Num_Mice = n_distinct(MouseID), .groups = "drop")
  
  # Render as table
  output$spont_fitc_4_summary_table <- renderTable({
    summary_fitc_4
  })
  
  
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
  
  
  icpms <- read.csv(here("data/phenotype/SLC_ICP-MS/Analysis_ICP_MS.csv"))
  icpms_long <- icpms %>% 
    filter(SampleType=="MUC-SI") %>%
    pivot_longer(cols = 2:8, names_to = "Element", values_to = "Concentration")
  
  icpms <- make_icp_ms_plot(icpms_long, "ICP-MS Cohort", stat_comparisons = list(c("WT","MUT")))
  
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
  output$plot_slc_spont_1A <- renderPlot({ print(plot_slc_spont_1A) })
  output$plot_slc_spont_1B <- renderPlot({ print(plot_slc_spont_1B) })
  output$plot_slc_spont_2 <- renderPlot({ print(plot_slc_spont_2 ) })
  output$plot_slc_spont_3 <- renderPlot({ print(plot_slc_spont_3 ) })
  output$plot_slc_spont_1C <- renderPlot({ print(plot_slc_spont_1C) })
  output$plot_slc_spont_1D <- renderPlot({ print(plot_slc_spont_1D) })
  
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
  
  
  combined_matrix <- read.csv(here("data/proteomics/combined_matrix.csv"), row.names=1)
  names(combined_matrix) <- gsub("X","", names(combined_matrix))
  s_meta <- read.csv(here("data/proteomics/Proteomics_Data - Metadata.csv"))
  
  selenoproteins_only <- combined_matrix %>% 
    filter(str_detect(Protein.names,"seleno")| 
             str_detect(Protein.names,"glutathione") |
             str_detect(Protein.names,"thioredoxin") |
             str_detect(Gene.Names,"Slc39") |
             str_detect(Gene.Names, "Iodothyronine")) %>%
    pivot_longer(2:37, names_to = "SampleID", values_to = "intensity") %>% 
    left_join(s_meta, by= "SampleID") 
  
  heatmap_mat <- selenoproteins_only %>%
    #filter(Sex=="Male") %>%
    #filter(UniprotID %in% top_proteins) %>%
    select(UniprotID, SampleID, intensity) %>%
    pivot_wider(names_from = SampleID, values_from = intensity) %>%
    column_to_rownames("UniprotID") %>%
    as.matrix()
  
  annotation_col <- s_meta %>%
    select(Genotype, Sex,Batch) %>%
    data.frame(row.names = s_meta$SampleID)
  
  annotation_row <- selenoproteins_only %>%
    dplyr::select(c("UniprotID", "Protein.names")) %>%
    unique() %>%
    column_to_rownames("UniprotID")
  
  heatmap <- pheatmap(
    heatmap_mat,
    annotation_col = annotation_col,
    #annotation_row= annotation_row,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
  )
  
  output$plot_slc_spont_proteomics <- renderPlot({
    print(heatmap)
  })
  # output$plot_slc_spont_proteomics <- renderPlot({ print(heatmap)})
  
  output$plot_smt_indo_1 <- renderPlot({ print(plot_smt_indo_1) })
  output$plot_slc_indo_1 <- renderPlot({ print(plot_slc_indo_1) })
  output$plot_slc_indo_2 <- renderPlot({ print(plot_slc_indo_2) })
  
  
  
  output$plot_icp_1 <- renderPlot({ print(icpms) })
  
  # Use the helper inside your outputs
  output$gsea_manhattan_mut <- renderPlotly({
    make_gsea_plot("results/RNA_seq/GSEA/M2_GSEA_HFD_Positive_MUT_vs_WT.csv")
  })
  
  output$gsea_manhattan_het <- renderPlotly({
    make_gsea_plot("results/RNA_seq/GSEA/M2_GSEA_HFD_Positive_HET_vs_WT.csv")
  })
  

  # Source modular server files
  files_to_source <- list.files(here("Server_Code"), pattern = "\\.R$", ignore.case = TRUE)
  
  cat("\n\n IN SERVER.R. Sourcing ... \n\n")
  
  for (file in files_to_source) {
    print(file)
    source(
      paste0("Server_Code/", file),
      local = TRUE
    )
  }
  
  
}
