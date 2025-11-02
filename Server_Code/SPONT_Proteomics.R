

observeEvent(input$proteomics_folder, {  # watch the selector
  req(input$proteomics_folder)           # ensure a value is selected
    
    # Map selection to actual folder paths
    data_dir <- switch(input$proteomics_folder,
                       "2025_only" = here("results/proteomics/2025_only/"),
                       "Combined"  = here("results/proteomics/combined/"))
    
    # List CSV files
    csv_files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
    if(length(csv_files) == 0) {
      output$tables_ui <- renderUI(h6("No CSV files found in this folder."))
      return()
    }
    
    
    # Read them all into a named list
    data_list <- lapply(csv_files, read.csv, stringsAsFactors = FALSE)
    names(data_list) <- basename(csv_files)
    
    # Dynamically create a datatable for each
    output$tables_ui <- renderUI({
      table_outputs <- lapply(names(data_list), function(name) {
        tagList(
          h6(name),
          DTOutput(outputId = paste0("table_", name)),
          hr()
        )
      })
      do.call(tagList, table_outputs)
    })
    
    # Render each datatable
    lapply(names(data_list), function(name) {
      output[[paste0("table_", name)]] <- renderDT({
        datatable(
          data_list[[name]],
          class = "display nowrap compact truncate-table",  # keeps visual style + local class
          escape = FALSE,
          options = list(
            pageLength = 5,
            scrollX = TRUE,
            # columnDefs: apply a JS renderer to all columns to produce truncated HTML with a title
            columnDefs = list(
              list(
                targets = "_all",
                render = JS(
                  "function(data, type, row, meta) {",
                  "  if(type !== 'display' || data === null) return data;",
                  "  var text = String(data);",
                  "  var safe = $('<div/>').text(text).html();",           
                  "  var max = 50; /* chars to show before truncation - change if desired */",
                  "  if(safe.length > max) {",
                  "    var short = safe.substring(0, max) + '\\u2026';", 
                  "    return '<div class=\"cell-truncate\" title=\"' + safe + '\">' + short + '</div>';",
                  "  }",
                  "  return '<div class=\"cell-truncate\" title=\"' + safe + '\">' + safe + '</div>';",
                  "}"
                )
              )
            )
          )
        )
      })
    })
    
  }
)