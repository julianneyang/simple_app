# SMT Negative Tab

observeEvent(input$tabs, {
  if (input$tabs == "SMT_Negative_Metagenomics") {
    species <- read.csv(here("data/sequencing_data/Metagenomics/Filtered_SMT_Shotgun_counts.csv")) %>% 
      group_by(X) %>%
      summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
      column_to_rownames("X")
    
    meta <-read.csv(here("data/sequencing_data/Metagenomics/Filtered_SMT_Shotgun_metadata.csv"), row.names=1) %>% 
      rownames_to_column("SampleID")
    
    
    ### Stacked Column Chart ---
    ld <- species %>% 
      rownames_to_column("Species") %>%
      pivot_longer(cols = starts_with("SLC"),names_to = "SampleID", values_to = "Count")
    
    ld <- ld %>% 
      group_by(SampleID) %>%
      mutate(Percentage = Count / sum(Count) * 100)
    
    ld <- merge(meta, ld, by= "SampleID")
    
    ld <- ld %>% 
      group_by(SampleID) %>%
      filter(Percentage >= 0.01) %>%
      mutate(Percentage = Count / sum(Count) * 100)
    
    ld$Genotype <- factor(ld$Genotype, levels=c("WT", "MUT"))
    
    cols_large <- unlist(lapply(1:nrow(palettes_d_names), function(i) {
      package_name <- palettes_d_names$package[i]
      palette_name <- palettes_d_names$palette[i]
      num_colors <- palettes_d_names$length[i]
      
      palette_colors <- paletteer_d(paste0(package_name, "::", palette_name), num_colors)
      
      return(palette_colors)
    }))
    set.seed(11)
    scrambled_cols <- unique(sample(cols_large))
    
    p <- ggplot(
      ld,
      aes(
        x = SampleID,
        y = Percentage,
        fill = Species,
        text = paste(
          "Species:", Species,
          "<br>Sample:", SampleID,
          "<br>Percent:", round(Percentage, 2)
        )
      )
    ) +
      geom_bar(stat = "identity", position = "stack") +
      facet_wrap(~Genotype, scales = "free_x") +
      labs(
        x = "Sample",
        y = "Species (%)",
        title = "Relative Abundance"
      ) +
      theme_cowplot(12) +
      scale_fill_manual(name = "Species", values = scrambled_cols) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    output$smt_microbiome_barplot  <- renderPlotly({
      ggplotly(p, tooltip = "text")
    }) 
    
    ### Beta Diversity ---
    
    source(here("src/Functions.R"))

    0.15*50
    filt_species <- prevalence_filter(species,8)
    filt_species <- filt_species[, colSums(filt_species) != 8]

    species.dist <- calculate_rsjensen(filt_species)

    cols <- c("WT"="black", "HET"="navy", "MUT"="firebrick")
    smt_pcoa <- generate_pcoA_plots(distance_matrix=species.dist,
                                counts = filt_species,
                                metadata = meta,
                                title="SMT_species",
                                colorvariable = Genotype,
                                colorvector = cols,
                                wa_scores_filepath = here("results/SMT_Maaslin2/pca_scores.csv"))
    output$smt_pcoa_plot <- renderPlot({
      smt_pcoa
    })
    
    output$labeled_smt_pcoa_plot <- renderPlot({
      smt_pcoa + aes(label=MouseID) + geom_label()
    })
    
    ### Plotting the significant Species ---
    
    significant_results <- read.delim(here("results/SMT_Maaslin3/Full/Sequencing_Depth_Site_Sex_Genotype/significant_results.tsv")) %>%
      filter(metadata=="Genotype") %>%
      filter(coef!="NA")
    
    abundance <- significant_results %>% 
      filter(null_hypothesis!=0)
    
    prevalence <- significant_results %>% 
      filter(null_hypothesis==0)
    
    make_taxa_dotplot <- function(ASV_significant_results_dataset,
                                  titlestring, 
                                  colorvector=  c("WT"="black", "MUT" = "firebrick"), qvalue=0.05){
      #data <- lc_dat_mut
      data <- as.data.frame(ASV_significant_results_dataset)
      data$annotation <- data$feature
      res_plot <- data %>% select(c("coef", "qval_individual","annotation"))
      res_plot <- unique(res_plot)
      res_plot <- res_plot %>%
        mutate(site = ifelse(coef< 0, "WT", "MUT"))
      
      y = tapply(res_plot$coef, res_plot$annotation, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
      y = sort(y, FALSE)   #switch to TRUE to reverse direction
      res_plot$annotation= factor(as.character(res_plot$annotation), levels = names(y))
      
      
      g1 <- res_plot %>%
        arrange(coef) %>%
        filter(qval_individual < qvalue, abs(coef) > 0) %>%
        ggplot(aes(
          x = coef,
          y = annotation,
          fill = site,
          text = paste0(
            "Effect size: ", round(coef, 3), "<br>",
            "q-value: ", signif(qval_individual, 3)
          )
        )) +
        geom_bar(stat = "identity") +
        cowplot::theme_cowplot(12) +
        theme(axis.text.y = element_text(face = "italic")) +
        scale_fill_manual(values = colorvector) +
        labs(x = "Effect size (MUT/WT)", y = "", fill = "") +
        theme(legend.position = "none") +
        ggtitle({{titlestring}}) +
        theme(plot.title = element_text(hjust = 0.5))
      
      # Convert to plotly with hover tooltips
      g1_plotly <- ggplotly(g1, tooltip = "text")
      
      return(g1_plotly)
      
    }
    
   smt_abundance_plot <- make_taxa_dotplot(ASV_significant_results_dataset = abundance,
                      #Relative_Abundance_filepath_rds = "results/SMT_Maaslin2/Relative_Abundance_Species_Luminal_SI_ASV.RDS",
                      titlestring = "Abundance: SMT Ile + Jej (MUT vs WT) ",
                      qvalue = 0.25)
   smt_prevalence_plot <- make_taxa_dotplot(ASV_significant_results_dataset = prevalence,
                      #Relative_Abundance_filepath_rds = "results/SMT_Maaslin2/Relative_Abundance_Species_Luminal_SI_ASV.RDS",
                      titlestring = "Prevalence: SMT Ile + Jej (MUT vs WT) ",
                      qvalue = 0.25)
   
   output$smt_DAT <- renderPlotly({
     smt_abundance_plot
   })
    
   output$smt_DPT <- renderPlotly({
     smt_prevalence_plot
   })
   
  }
})