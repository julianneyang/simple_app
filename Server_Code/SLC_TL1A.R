# SLC TL1A Tab

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