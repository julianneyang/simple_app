# SMT Negative Tab

observeEvent(input$tabs, {
  if (input$tabs == "SMT_Negative_Metagenomics") {
    species <- read.csv(here("data/sequencing_data/Metagenomics/Filtered_SMT_Shotgun_counts.csv")) %>% 
      group_by(X) %>%
      summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
      column_to_rownames("X")
    
    meta <-read.csv(here("data/sequencing_data/Metagenomics/Filtered_SMT_Shotgun_metadata.csv"), row.names=1) %>% 
      rownames_to_column("SampleID")
    
    
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
    
    
  }
})