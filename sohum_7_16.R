# Needed libraries for the application
# install.packages(c("shiny", "dplyr", "tidyr", "readr", "ggplot2", "plotly", "patchwork"))

library(shiny)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(plotly) # Added for interactivity
library(patchwork)

# Set Shiny's upload limit to 600 MB
options(shiny.maxRequestSize = 600 * 1024^2)

# --- User Interface (UI) ---
ui <- fluidPage(
  titlePanel("Interactive Gene Expression Visualization Tool"),
  sidebarLayout(
    sidebarPanel(
      fileInput("expr", "1. Expression Matrix (.tsv or .tsv.gz)",
                accept = c(".tsv", ".tsv.gz")),
      fileInput("meta", "2. Metadata File (.tsv or .tsv.gz)",
                accept = c(".tsv", ".tsv.gz")),
      hr(),
      uiOutput("gene_ui"),
      uiOutput("group_ui"),
      uiOutput("split_ui"),
      p(strong("Note:"), "Select at least one gene to generate plots.")
    ),
    mainPanel(
      tabsetPanel(
        # Changed plotOutput to plotlyOutput for all tabs
        tabPanel("Feature Plot", plotlyOutput("featPlot", height = "600px")),
        tabPanel("Bubble Plot", plotlyOutput("dotPlot", height = "800px")),
        tabPanel("Violin Plot", plotlyOutput("vlnPlot", height = "600px")),
        tabPanel("Box Plot", plotlyOutput("boxPlot", height = "600px")),
        tabPanel("Heatmap", plotlyOutput("heatmapPlot", height = "800px"))
      )
    )
  )
)

# --- Server Logic ---
server <- function(input, output, session) {
  
  # Reactive expression to load and merge data
  plot_data <- reactive({
    req(input$expr, input$meta)
    meta_df <- read_tsv(input$meta$datapath) %>%
      tibble::column_to_rownames(var = colnames(.)[1])
    expr_df <- read_tsv(input$expr$datapath)
    genes <- expr_df[[1]]
    mat <- as.matrix(expr_df[,-1])
    rownames(mat) <- genes
    t_mat <- t(mat)
    ordered_meta <- meta_df[rownames(t_mat), ]
    # Add cell names as a column for easier use in ggplot
    combined_df <- as.data.frame(t_mat) %>%
      mutate(cell = rownames(.)) %>%
      bind_cols(ordered_meta)
    return(combined_df)
  })
  
  # Dynamic UI rendering (no changes here)
  output$gene_ui <- renderUI({
    req(plot_data())
    all_cols <- colnames(plot_data())
    meta_cols <- c("cell", colnames(read_tsv(input$meta$datapath) %>% select(-1)))
    gene_choices <- setdiff(all_cols, meta_cols)
    selectInput("genes", "Select Gene(s):", choices = gene_choices, selected = gene_choices[1], multiple = TRUE)
  })
  output$group_ui <- renderUI({
    req(plot_data())
    meta_cols <- colnames(read_tsv(input$meta$datapath) %>% select(-1))
    group_choices <- grep("UMAP", meta_cols, value = TRUE, invert = TRUE)
    selectInput("group", "Group By (X-axis/Categories):", choices = group_choices, selected = group_choices[1])
  })
  output$split_ui <- renderUI({
    req(plot_data())
    meta_cols <- colnames(read_tsv(input$meta$datapath) %>% select(-1))
    split_choices <- c("None", grep("UMAP", meta_cols, value = TRUE, invert = TRUE))
    selectInput("split", "Split By (Rows/Panels):", choices = split_choices, selected = "None")
  })
  
  # --- Interactive Plotting Functions ---
  
  # 1. Feature Plot (UMAP)
  output$featPlot <- renderPlotly({
    req(input$genes, plot_data(), "UMAP_Xaxis" %in% names(plot_data()))
    df <- plot_data()
    plots <- lapply(input$genes, function(gene) {
      p <- ggplot(df, aes(x = UMAP_Xaxis, y = UMAP_Yaxis, color = .data[[gene]], text = paste0("Cell: ", cell, "<br>Expr: ", round(.data[[gene]], 2)))) +
        geom_point(size = 0.5, alpha = 0.7) +
        labs(title = gene, color = "Expression") +
        scale_color_viridis_c() +
        theme_minimal()
      if (input$split != "None") {
        p <- p + facet_wrap(vars(.data[[input$split]]))
      }
      p
    })
    ggplotly(wrap_plots(plots, ncol = 1), tooltip = "text")
  })
  
  # 2. Bubble Plot
  output$dotPlot <- renderPlotly({
    req(input$genes, input$group, plot_data())
    dot_plot_data <- plot_data() %>%
      select(all_of(c(input$genes, input$group, if(input$split != "None") input$split))) %>%
      pivot_longer(cols = all_of(input$genes), names_to = "gene", values_to = "expression")
    summary_df <- dot_plot_data %>%
      group_by(gene, .data[[input$group]], .add = TRUE) %>%
      { if (input$split != "None") group_by(., .data[[input$split]], .add = TRUE) else . } %>%
      summarise(avg_expr = mean(expression), pct_expr = (sum(expression > 0) / n()) * 100) %>%
      ungroup()
    
    p <- ggplot(summary_df, aes(x = gene, y = .data[[input$group]], size = pct_expr, color = avg_expr,
                                text = paste0("Avg. Expr: ", round(avg_expr, 2), "<br>% Expr: ", round(pct_expr, 2), "%"))) +
      geom_point(alpha = 0.8) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Avg. Expr.") +
      scale_size_continuous(name = "% Expr.") +
      labs(x = "Feature", y = "Cell Subset") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    if (input$split != "None") {
      p <- p + facet_grid(rows = vars(.data[[input$split]]), scales = "free_y", space = "free_y")
    }
    ggplotly(p, tooltip = "text")
  })
  
  # 3. Violin Plot
  output$vlnPlot <- renderPlotly({
    req(input$genes, input$group, plot_data())
    
    p <- ggplot(plot_data(), aes(x = .data[[input$group]], y = .data[[input$genes[1]]], fill = .data[[input$group]])) +
      # Add jitter points for individual cells
      geom_jitter(size = 0.4, alpha = 0.5, height = 0, width = 0.2) +
      # Add the violin plot layer
      geom_violin(trim = FALSE, scale = "width", alpha = 0.7) +
      # Add the boxplot inside the violin
      geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.7) +
      labs(title = paste("Expression of", input$genes[1]), x = input$group, y = "Expression") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    
    if (input$split != "None") {
      p <- p + facet_wrap(vars(.data[[input$split]]))
    }
    ggplotly(p)
  })
  
  # 4. Box Plot
  output$boxPlot <- renderPlotly({
    req(input$genes, input$group, plot_data())
    p <- ggplot(plot_data(), aes(x = .data[[input$group]], y = .data[[input$genes[1]]], fill = .data[[input$group]])) +
      geom_boxplot() +
      labs(title = paste("Expression of", input$genes[1]), x = input$group, y = "Expression") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    if (input$split != "None") { p <- p + facet_wrap(vars(.data[[input$split]])) }
    ggplotly(p)
  })
  
  # 5. Heatmap (rewritten with ggplot for interactivity)
  output$heatmapPlot <- renderPlotly({
    req(length(input$genes) >= 2, input$group, plot_data())
    
    heatmap_data <- plot_data() %>%
      arrange(.data[[input$group]]) %>%
      mutate(cell = factor(cell, levels = unique(cell))) %>%
      select(cell, all_of(input$genes), all_of(input$group)) %>%
      pivot_longer(cols = all_of(input$genes), names_to = "gene", values_to = "expression") %>%
      group_by(gene) %>%
      mutate(scaled_expression = as.numeric(scale(expression))) %>%
      ungroup()
    
    p <- ggplot(heatmap_data, aes(x = cell, y = gene, fill = scaled_expression,
                                  text = paste0("Cell: ", cell, "<br>Gene: ", gene, "<br>Group: ", .data[[input$group]], "<br>Scaled Expr: ", round(scaled_expression, 2)))) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", name = "Scaled Expr.") +
      labs(x = "Cells", y = "Genes") +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()
      )
    
    ggplotly(p, tooltip = "text")
  })
}

# Launch the Shiny app
shinyApp(ui, server)
