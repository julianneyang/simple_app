# Selenium Supplementation Tab

observeEvent(input$tabs, {
  if (input$tabs == "SE Supplementation") {
    
    se_indo <- read.csv(here("data/phenotype/SE_Supp/SE_Supp_Indomethacin_DAI.csv"))
    se_indo_ctrl <- se_indo %>% filter(Diet =="Control")
    plot_se_supp_1 <- plot_avg_trajectory(se_indo_ctrl, title_string = "Indomethacin", subtitle_string = "Control Diet; HET p = 0.5727, MUT p = 0.5037")
    
    se_indo_supp <- se_indo %>% filter(Diet =="Selenium")
    plot_se_supp_2 <- plot_avg_trajectory(se_indo_supp, title_string = "Indomethacin", subtitle_string = "Selenium Diet; HET p = 0.7391, MUT p = 0.6525")
    
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
    
    se_supp_icpms <- read.csv(here("data/phenotype/SE_Supp/Se_Supp_ICP_MS_Analysis.csv"))
    se_supp_icpms_long <- se_supp_icpms %>% 
      filter(Weight>=3) %>%
      pivot_longer(cols = 3:9, names_to = "Element", values_to = "Concentration")
    control_diet <- se_supp_icpms_long %>% filter(Diet=="Control")
    se_supp_diet <- se_supp_icpms_long %>% filter(Diet == "Selenium")
    
    control_icpms <- make_icp_ms_plot(control_diet, "Control Diet", stat_comparisons = full_comparisons)
    
    se_supp_icpms <- make_icp_ms_plot(se_supp_diet,  "Se Supp Diet", stat_comparisons = full_comparisons)
    
    diet_icpms <- make_icp_ms_plot_diet_comparison(se_supp_icpms_long, "Within Diet")
    
    se_supp_fitc <- plot_fitc(here("data/phenotype/SE_Supp/Se_Supp_FITC_Running_Total_2025Oct.csv"), "FITC", "Selenium Diet; HET p= 0.5557, MUT p=0.0656") 
    
    control_fitc <- plot_fitc(here("data/phenotype/SE_Supp/Control_FITC_Running_Total_2025Oct.csv"), "FITC", "Control Diet; HET p = 0.3787, MUT p = 0.9314") 
    
    output$plot_se_supp_1 <- renderPlot({ print(plot_se_supp_1) })
    output$plot_se_supp_2 <- renderPlot({ print(plot_se_supp_2) })
    output$plot_se_supp_3 <- renderPlot({ print(se_histo_ctrl) })
    output$plot_se_supp_4 <- renderPlot({ print(se_histo_indo) })
    output$plot_se_supp_5 <- renderPlot({ print(control_icpms) })
    output$plot_se_supp_6 <- renderPlot({ print(se_supp_icpms) })
    output$plot_se_supp_7 <- renderPlot({ print(control_fitc) })
    output$plot_se_supp_8 <- renderPlot({ print(se_supp_fitc) })
    output$plot_se_supp_9 <- renderPlot({ print(diet_icpms) })
  }
})
