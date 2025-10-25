# SMT Negative Tab

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