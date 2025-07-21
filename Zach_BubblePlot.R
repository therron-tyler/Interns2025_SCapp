# Cleaned and accessible version of your Shiny gene expression app
# Functionalities remain the same, with improvements in accessibility and UI

library(shiny)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(plotly)
library(patchwork)

options(shiny.maxRequestSize = 600 * 1024^2)
Sys.setenv("VROOM_CONNECTION_SIZE" = 10000000)

ui <- fluidPage(
  titlePanel("Interactive Gene Expression Visualization Tool"),
  sidebarLayout(
    sidebarPanel(
      fileInput("expr", "1. Expression Matrix (.tsv or .tsv.gz)", accept = c(".tsv", ".tsv.gz")),
      fileInput("meta", "2. Metadata File (.tsv or .tsv.gz)", accept = c(".tsv", ".tsv.gz")),
      hr(),
      uiOutput("gene_ui"),
      uiOutput("group_ui"),
      uiOutput("split_ui"),
      p(strong("Note:"), "Select at least one gene to generate plots."),
      tags$style("label { font-weight: bold; }")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Feature Plot", plotlyOutput("featPlot", height = "600px")),
        tabPanel("Bubble Plot", plotlyOutput("dotPlot", height = "800px")),
        tabPanel("Violin Plot", plotlyOutput("vlnPlot", height = "600px")),
        tabPanel("Box Plot", plotlyOutput("boxPlot", height = "600px")),
        tabPanel("Heatmap", plotlyOutput("heatmapPlot", height = "800px"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  plot_data <- reactive({
    req(input$expr, input$meta)
    meta_df <- read_tsv(input$meta$datapath, show_col_types = FALSE) %>%
      tibble::column_to_rownames(var = colnames(.)[1])
    expr_df <- read_tsv(input$expr$datapath, show_col_types = FALSE)
    genes <- expr_df[[1]]
    mat <- as.matrix(expr_df[,-1])
    rownames(mat) <- genes
    t_mat <- t(mat)
    ordered_meta <- meta_df[rownames(t_mat), , drop = FALSE]
    
    df <- as.data.frame(t_mat) %>%
      mutate(cell = rownames(.)) %>%
      bind_cols(ordered_meta)
    
    df
  })
  
  output$gene_ui <- renderUI({
    req(plot_data())
    all_cols <- colnames(plot_data())
    meta_cols <- c("cell", colnames(read_tsv(input$meta$datapath, show_col_types = FALSE)[-1]))
    gene_choices <- setdiff(all_cols, meta_cols)
    selectInput("genes", "Select Gene(s):", choices = gene_choices, selected = gene_choices[1], multiple = TRUE)
  })
  
  output$group_ui <- renderUI({
    req(plot_data())
    meta_cols <- colnames(read_tsv(input$meta$datapath, show_col_types = FALSE)[-1])
    group_choices <- meta_cols
    selectInput("group", "Group By (X-axis/Categories):", choices = group_choices, selected = group_choices[1])
  })
  
  output$split_ui <- renderUI({
    req(plot_data())
    meta_cols <- colnames(read_tsv(input$meta$datapath, show_col_types = FALSE)[-1])
    split_choices <- c("None", meta_cols)
    selectInput("split", "Split By (Rows/Panels):", choices = split_choices, selected = "None")
  })
  
  output$featPlot <- renderPlotly({
    req(input$genes, plot_data())
    df <- plot_data()
    umap_cols <- grep("UMAP", colnames(df), value = TRUE)
    req(length(umap_cols) >= 2)
    
    xcol <- umap_cols[1]
    ycol <- umap_cols[2]
    
    plots <- lapply(input$genes, function(gene) {
      p <- ggplot(df, aes(x = .data[[xcol]], y = .data[[ycol]], color = .data[[gene]],
                          text = paste0("Cell: ", cell, "<br>Expr: ", round(.data[[gene]], 2)))) +
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
  
  output$dotPlot <- renderPlotly({
    req(input$genes, input$group, plot_data())
    df <- plot_data()
    
    long_df <- df %>%
      select(all_of(c("cell", input$genes, input$group, if (input$split != "None") input$split))) %>%
      pivot_longer(cols = all_of(input$genes), names_to = "gene", values_to = "expression")
    
    summary_df <- long_df %>%
      group_by(gene, .data[[input$group]]) %>%
      { if (input$split != "None") group_by(., .data[[input$split]], .add = TRUE) else . } %>%
      summarise(avg_expr = mean(expression), pct_expr = (sum(expression > 0) / n()) * 100, .groups = "drop")
    
    p <- ggplot(summary_df, aes(x = gene, y = .data[[input$group]], size = pct_expr, color = avg_expr,
                                text = paste0("Avg Expr: ", round(avg_expr, 2), "<br>% Expr: ", round(pct_expr, 1), "%"))) +
      geom_point(alpha = 0.8) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      scale_size_continuous(name = "% Expr.") +
      labs(x = "Gene", y = "Group") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    if (input$split != "None") {
      p <- p + facet_grid(rows = vars(.data[[input$split]]), scales = "free_y", space = "free_y")
    }
    
    ggplotly(p, tooltip = "text")
  })
  
  output$vlnPlot <- renderPlotly({
    req(input$genes, input$group, plot_data())
    df <- plot_data()
    gene <- input$genes[1]
    
    p <- ggplot(df, aes(x = .data[[input$group]], y = .data[[gene]], fill = .data[[input$group]])) +
      geom_violin(trim = FALSE, scale = "width", alpha = 0.7) +
      geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.7) +
      geom_jitter(size = 0.4, alpha = 0.5, width = 0.2) +
      labs(title = paste("Expression of", gene), y = "Expression") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    
    if (input$split != "None") {
      p <- p + facet_wrap(vars(.data[[input$split]]))
    }
    
    ggplotly(p)
  })
  
  output$boxPlot <- renderPlotly({
    req(input$genes, input$group, plot_data())
    df <- plot_data()
    gene <- input$genes[1]
    
    p <- ggplot(df, aes(x = .data[[input$group]], y = .data[[gene]], fill = .data[[input$group]])) +
      geom_boxplot() +
      labs(title = paste("Expression of", gene), y = "Expression") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    
    if (input$split != "None") {
      p <- p + facet_wrap(vars(.data[[input$split]]))
    }
    
    ggplotly(p)
  })
  
  output$heatmapPlot <- renderPlotly({
    req(length(input$genes) >= 2, input$group, plot_data())
    df <- plot_data()
    
    heatmap_data <- df %>%
      arrange(.data[[input$group]]) %>%
      mutate(cell = factor(cell, levels = unique(cell))) %>%
      select(cell, all_of(input$genes), all_of(input$group)) %>%
      pivot_longer(cols = all_of(input$genes), names_to = "gene", values_to = "expression") %>%
      group_by(gene) %>%
      mutate(scaled_expression = as.numeric(scale(expression))) %>%
      ungroup()
    
    p <- ggplot(heatmap_data, aes(x = cell, y = gene, fill = scaled_expression,
                                  text = paste0("Cell: ", cell, "<br>Gene: ", gene,
                                                "<br>Group: ", .data[[input$group]],
                                                "<br>Scaled Expr: ", round(scaled_expression, 2)))) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", name = "Scaled Expr.") +
      labs(x = "Cells", y = "Genes") +
      theme_minimal() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())
    
    ggplotly(p, tooltip = "text")
  })
}

shinyApp(ui, server)
