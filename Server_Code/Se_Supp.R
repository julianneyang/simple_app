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
    
    se_fitc_df <- read.csv(here("data/phenotype/SE_Supp/Se_Supp_FITC_Running_Total_2025Oct.csv")) %>%
      filter(Omit=="No")
    se_supp_fitc <- plot_fitc(dataframe = se_fitc_df, "FITC", "Selenium Diet; HET p= 0.5557, MUT p=0.0656") 
    
    control_fitc_df <- read.csv(here("data/phenotype/SE_Supp/Control_FITC_Running_Total_2025Oct.csv")) %>% 
      filter(Omit =="No")
    control_fitc <- plot_fitc(dataframe=control_fitc_df, "FITC", "Control Diet; HET p = 0.3787, MUT p = 0.9314") 
    
    merged_fitc <- rbind(se_fitc_df, control_fitc_df) %>%
      mutate(Diet = str_replace_all(Diet, "Selenium", "Se")) %>%
      mutate(Diet = str_replace_all(Diet, "Control", "C")) %>%
      mutate(Diet_Genotype = paste0(Diet,"_",Genotype))
    
    merged_fitc$Diet_Genotype <- factor(merged_fitc$Diet_Genotype, levels=c(
      "C_WT", "Se_WT",
      "C_HET","Se_HET",
      "C_MUT", "Se_MUT"))
    
    merged_fitc_plot <- plot_fitc_diet_comparison(dataframe = merged_fitc, 
                             title_string = "FITC",
                             subtitle_string = "Se vs Ctrl Diet",
                             stat_comparisons = list( c("C_WT", "Se_WT"), 
                                                      c("C_HET","Se_HET"),
                                                      c("C_MUT", "Se_MUT")))
    
    control_spont_fitc <- control_fitc_df %>% 
      dplyr::select(c(MouseID, Genotype, Sex, Plasma_FITC,Size, Diet)) 
    
    spont_fitc_agg <- spont_fitc_df %>% 
      mutate(Diet="Control") %>% 
      filter(Omit=="No") %>% 
      dplyr::select(-c("Batch", "Omit"))
    
    merged_control_agg <- rbind(control_spont_fitc, spont_fitc_agg) 
    
    control_agg_plot <- plot_fitc(dataframe=merged_control_agg, "FITC", "CD and Spont; HET p = 0.044, MUT p = 0.000224") 
    
    se_fitc_df <- se_fitc_df %>%  dplyr::select(c(MouseID, Genotype, Sex, Plasma_FITC,Size, Diet)) 
    merged_full <- rbind(merged_control_agg, se_fitc_df )%>%
      mutate(Diet = str_replace_all(Diet, "Selenium", "Se")) %>%
      mutate(Diet = str_replace_all(Diet, "Control", "C")) %>%
      mutate(Diet_Genotype = paste0(Diet,"_",Genotype)) 
    merged_full$Diet_Genotype <- factor(merged_full$Diet_Genotype, levels=c(
      "C_WT", "Se_WT",
      "C_HET","Se_HET",
      "C_MUT", "Se_MUT"))
    
      
      
    
    merged_fitc_spont_plot <- plot_fitc_diet_comparison(dataframe = merged_full, 
                                                        title_string = "FITC Control and Spont",
                                                        subtitle_string = "Se vs Ctrl Diet",
                                                        stat_comparisons = list( c("C_WT", "Se_WT"), 
                                                                                 c("C_HET","Se_HET"),
                                                                                 c("C_MUT", "Se_MUT")))
      
    output$plot_se_supp_1 <- renderPlot({ print(plot_se_supp_1) })
    output$plot_se_supp_2 <- renderPlot({ print(plot_se_supp_2) })
    output$plot_se_supp_3 <- renderPlot({ print(se_histo_ctrl) })
    output$plot_se_supp_4 <- renderPlot({ print(se_histo_indo) })
    output$plot_se_supp_5 <- renderPlot({ print(control_icpms) })
    output$plot_se_supp_6 <- renderPlot({ print(se_supp_icpms) })
    output$plot_se_supp_7 <- renderPlot({ print(control_fitc) })
    output$plot_se_supp_8 <- renderPlot({ print(se_supp_fitc) })
    output$plot_se_supp_10 <- renderPlot({ print(merged_fitc_plot) })
    output$plot_se_supp_9 <- renderPlot({ print(diet_icpms) })
    
    output$plot_se_supp_11 <- renderPlot({ print(control_agg_plot) })
    output$plot_se_supp_12 <- renderPlot({ print(merged_fitc_spont_plot) })
  }
})
