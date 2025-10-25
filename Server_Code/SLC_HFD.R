# SLC HFD Tab 

observeEvent(input$tabs, {
  if (input$tabs == "SLC_HFD") {
    output$plot_slc_hfd_1 <- renderPlot({ print(plot_hfd_1) })
    output$plot_slc_hfd_2 <- renderPlot({ print(plot_hfd_2) })
    plot_hfd_1 <- plot_histology(csv_filepath= here("data/phenotype/SLC_HFD/HFD_Histo.csv"),
                                 title_string = "Ileum Histology",
                                 subtitle_string =  "HFD Cohort",
                                 stat_comparisons =  full_comparisons) 
    
    hfd_fitc <- read.csv(here("data/phenotype/SLC_HFD/HFD_FITC_Analysis.csv"))
    plot_hfd_2 <- plot_fitc(dataframe = hfd_fitc,
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