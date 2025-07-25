# -----------------------------------------------------------------------------
# Interactive Gene Expression Visualization Tool
# Collaborative sections clearly marked for easy editing by each team member
# -----------------------------------------------------------------------------

# Required libraries
library(shiny)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(plotly)
library(patchwork)
library(RColorBrewer)

# Increase upload limit
options(shiny.maxRequestSize = 600 * 1024^2)
Sys.setenv("VROOM_CONNECTION_SIZE" = 10000000)

# --- UI DEFINITION -----------------------------------------------------------
ui <- fluidPage(
  titlePanel("Interactive Gene Expression Visualization Tool - Sohum, Suki, Kim, and Zach"),
  sidebarLayout(
    sidebarPanel(
      fileInput("expr", "1. Expression Matrix (.tsv or .tsv.gz)", accept = c(".tsv", ".tsv.gz")),
      fileInput("meta", "2. Metadata File (.tsv or .tsv.gz)",      accept = c(".tsv", ".tsv.gz")),
      hr(),
      uiOutput("gene_ui"),      # Gene selector
      uiOutput("group_ui"),     # Group-by selector
      uiOutput("split_ui"),     # Split-by selector
      checkboxInput("color_by_gene", 
                    "Feature Plot: Color by gene expression", 
                    value = FALSE),
      p(strong("Note:"), "For Violin/Box plots, multiple genes will be shown in a grid.")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Feature Plot", plotlyOutput("featPlot",  height = "600px")),
        tabPanel("Bubble Plot",  plotlyOutput("dotPlot",   height = "800px")),
        tabPanel("Violin Plot",  plotlyOutput("vlnPlot",   height = "800px")),
        tabPanel("Heatmap",      plotlyOutput("heatmapPlot",height = "800px")),
        tabPanel("Box Plot",     plotlyOutput("boxPlot",    height = "800px"))  # ← last
      )
    )
  )
)

# --- SERVER LOGIC -----------------------------------------------------------
server <- function(input, output, session) {
  
  # 1) Load and merge expression + metadata -------------------------------
  plot_data <- reactive({
    req(input$expr, input$meta)
    meta_df <- read_tsv(input$meta$datapath, show_col_types = FALSE) %>%
      tibble::column_to_rownames(var = colnames(.)[1])
    expr_df <- read_tsv(input$expr$datapath, show_col_types = FALSE)
    genes <- expr_df[[1]]
    mat   <- as.matrix(expr_df[,-1]); rownames(mat) <- genes
    df    <- t(mat) %>% as.data.frame() %>% tibble::rownames_to_column("cell")
    bind_cols(df, meta_df[df$cell, , drop = FALSE])
  })
  
  # 2) Dynamic UI for gene, group, split ----------------------------------
  output$gene_ui <- renderUI({
    req(plot_data())
    all_cols  <- colnames(plot_data())
    meta_cols <- c("cell", colnames(read_tsv(input$meta$datapath, show_col_types = FALSE)[-1]))
    gene_choices <- setdiff(all_cols, meta_cols)
    selectInput("genes", "Select Gene(s):",
                choices  = gene_choices,
                selected = gene_choices[1],
                multiple = TRUE)
  })
  output$group_ui <- renderUI({
    req(plot_data())
    meta_cols <- colnames(read_tsv(input$meta$datapath, show_col_types = FALSE)[-1])
    selectInput("group", "Group By (Categories/Colors):",
                choices  = meta_cols,
                selected = meta_cols[1])
  })
  output$split_ui <- renderUI({
    req(plot_data())
    meta_cols <- colnames(read_tsv(input$meta$datapath, show_col_types = FALSE)[-1])
    selectInput("split", "Split By (Rows/Panels):",
                choices  = c("None", meta_cols),
                selected = "None")
  })
  
  # === FEATURE PLOT SECTION START ===
  output$featPlot <- renderPlotly({
    req(plot_data(), "UMAP_Xaxis" %in% names(plot_data()), input$group)
    df <- plot_data()
    if (!input$color_by_gene) {
      p <- ggplot(df,
                  aes(x = UMAP_Xaxis, y = UMAP_Yaxis,
                      color = .data[[input$group]],
                      text = paste0("Cell: ", cell, "<br>Group: ", .data[[input$group]]))) +
        geom_point(size = 1, alpha = 0.8) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        labs(title = "UMAP Colored by Metadata Group", color = input$group) +
        theme_minimal()
    } else {
      req(input$genes)
      if (length(input$genes) == 1) {
        p <- ggplot(df,
                    aes(x = UMAP_Xaxis, y = UMAP_Yaxis,
                        color = .data[[input$genes]],
                        text  = paste0("Cell: ", cell, "<br>Expr: ", round(.data[[input$genes]],2)))) +
          geom_point(size = 1, alpha = 0.8) +
          scale_color_viridis_c() +
          labs(title = paste("Expression of", input$genes), color = "Expression") +
          theme_minimal()
      } else {
        df$ModuleScore <- rowMeans(df[, input$genes, drop = FALSE], na.rm = TRUE)
        p <- ggplot(df,
                    aes(x = UMAP_Xaxis, y = UMAP_Yaxis,
                        color = ModuleScore,
                        text  = paste0("Cell: ", cell, "<br>Score: ", round(ModuleScore,2)))) +
          geom_point(size = 1, alpha = 0.8) +
          scale_color_viridis_c() +
          labs(title = "Module Score for Selected Genes", color = "Avg. Expr.") +
          theme_minimal()
      }
    }
    ggplotly(p, tooltip = "text")
  })
  # === FEATURE PLOT SECTION END ===
  
  # === BUBBLE PLOT SECTION START ===
  output$dotPlot <- renderPlotly({
    req(input$genes, input$group, plot_data())
    df <- plot_data()
    
    # Use input$split as the group (y-axis)
    # Use input$group as the facet (split)
    split_var <- if (input$split != "None") input$split else NULL
    group_var <- input$group
    
    # Create long-format dataframe
    dot_data <- df %>%
      select(all_of(c(input$genes, split_var, group_var))) %>%
      pivot_longer(cols = all_of(input$genes), names_to = "gene", values_to = "expression")
    
    # Group and summarize
    summary_df <- dot_data %>%
      group_by(gene, .data[[split_var]], .add = TRUE) %>%
      { if (!is.null(group_var)) group_by(., .data[[group_var]], .add = TRUE) else . } %>%
      summarise(
        avg_expr = mean(expression, na.rm = TRUE),
        pct_expr = 100 * sum(expression > 0, na.rm = TRUE) / n(),
        .groups = "drop"
      )
    
    # Base plot
    p <- ggplot(summary_df, aes(
      x = gene,
      y = .data[[split_var]],
      size = pct_expr,
      color = avg_expr,
      text = paste0(
        "Gene: ", gene, "<br>",
        split_var, ": ", .data[[split_var]], "<br>",
        if (!is.null(group_var)) paste0(group_var, ": ", .data[[group_var]], "<br>") else "",
        "Avg. Expr: ", round(avg_expr, 2), "<br>",
        "% Expr: ", round(pct_expr, 1), "%"
      )
    )) +
      geom_point(alpha = 0.9) +
      geom_text(aes(label = paste0(round(pct_expr, 1), "%")),
                size = 2.7, vjust = -1.3, color = "black", show.legend = FALSE) +
      scale_color_viridis_c(name = "Avg. Expr.") +
      scale_size_continuous(name = "% Expr.", range = c(2, 10)) +
      labs(
        title = "Gene Expression Bubble Plot",
        x = "Gene",
        y = split_var
      ) +
      theme_minimal(base_size = 13) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5)
      )
    
    # Facet by the original "Group By" if applicable
    if (!is.null(group_var)) {
      p <- p + facet_wrap(vars(.data[[group_var]]), ncol = 1)
    }
    
    ggplotly(p, tooltip = "text") %>% layout(legend = list(orientation = "h", y = -0.2))
  })
  # === HEATMAP SECTION END ===
  
  # === BOX PLOT SECTION START ===
  output$boxPlot <- renderPlotly({
    req(input$genes, input$group, plot_data())
    df <- plot_data() %>%
      pivot_longer(cols = all_of(input$genes), names_to="gene", values_to="expr")
    p <- ggplot(df, aes(x=.data[[input$group]], y=expr, fill=.data[[input$group]])) +
      geom_boxplot() +
      facet_wrap(~gene, scales="free_y", ncol=1) +
      theme_bw() +
      theme(axis.text.x=element_text(angle=45,hjust=1),
            legend.position="none") +
      labs(x=input$group, y="Expression")
    if(input$split!="None") {
      p <- p + aes(fill=.data[[input$split]]) +
        theme(legend.position="top")
    }
    ggplotly(p)
  })
  # === BOX PLOT SECTION END ===
}

# Launch the app
shinyApp(ui, server)