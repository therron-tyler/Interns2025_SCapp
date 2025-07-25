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
    dot <- df %>%
      select(all_of(c(input$genes, input$group,
                      if(input$split!="None") input$split))) %>%
      pivot_longer(cols = all_of(input$genes), names_to="gene", values_to="expr")
    summary_df <- dot %>%
      group_by(gene, .data[[input$group]], .add=TRUE) %>%
      { if(input$split!="None") group_by(., .data[[input$split]], .add=TRUE) else . } %>%
      summarise(
        avg_expr = mean(expr),
        pct_expr = 100 * sum(expr>0)/n(),
        .groups="drop"
      )
    xcol <- if(input$split!="None") input$split else input$group
    ycol <- input$group
    p <- ggplot(summary_df,
                aes(x = .data[[xcol]], y = .data[[ycol]],
                    size = pct_expr, color = avg_expr)) +
      geom_point(alpha=0.8) +
      geom_text(aes(label=paste0(round(pct_expr,1),"%")),
                color="black", size=3, vjust=-1, show.legend=FALSE) +
      scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0) +
      scale_size_continuous(range=c(2,12)) +
      labs(x=xcol, y=ycol, title="Bubble Plot") +
      theme_minimal() +
      theme(axis.text.x=element_text(angle=45,hjust=1))
    p <- p + aes(text=paste0(
      "Gene: ", gene, "<br>",
      ycol, ": ", .data[[ycol]], "<br>",
      if(input$split!="None") paste0(xcol, ": ", .data[[xcol]], "<br>") else "",
      "Avg Expr: ", round(avg_expr,2), "<br>",
      "% Expr: ", round(pct_expr,2),"%"
    ))
    ggplotly(p, tooltip="text")
  })
  # === BUBBLE PLOT SECTION END ===
  
  # === VIOLIN PLOT SECTION START ===
  output$vlnPlot <- renderPlotly({
    req(input$genes, input$group, plot_data())
    df <- plot_data() %>%
      pivot_longer(cols = all_of(input$genes), names_to="gene", values_to="expr")
    p <- ggplot(df, aes(x=.data[[input$group]], y=expr)) +
      geom_violin(trim=FALSE, scale="width", alpha=0.7) +
      geom_boxplot(width=0.1, fill="white", outlier.shape=NA, alpha=0.5) +
      geom_jitter(size=0.4, alpha=0.3, width=0.2) +
      facet_wrap(~gene, scales="free_y", ncol=1) +
      theme_bw() +
      theme(axis.text.x=element_text(angle=45,hjust=1)) +
      labs(x=input$group,y="Expression")
    if(input$split!="None") {
      p <- p + aes(fill=.data[[input$split]]) +
        theme(legend.position="top")
    }
    ggplotly(p)
  })
  # === VIOLIN PLOT SECTION END ===
  
  # === HEATMAP SECTION START ===
  output$heatmapPlot <- renderPlotly({
    req(input$genes, input$group, plot_data())
    df <- plot_data()
    grp <- df[[input$group]]
    if(input$split!="None") grp <- paste(df[[input$split]],grp,sep="_")
    mat <- as.matrix(df[input$genes])
    avg_mat <- t(rowsum(mat, grp)/as.vector(table(grp)))
    mat_z <- t(scale(t(avg_mat))); mat_z[is.na(mat_z)]<-0
    plot_ly(
      x=colnames(mat_z), y=rownames(mat_z), z=mat_z,
      type="heatmap",
      colors=colorRamp(c("#440154","white","#21918c")),
      showscale=TRUE
    ) %>%
      layout(
        title="Expression Heatmap",
        xaxis=list(title="Group",tickangle=45),
        yaxis=list(title="Gene"),
        margin=list(l=80,r=20,b=80,t=50)
      )
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
