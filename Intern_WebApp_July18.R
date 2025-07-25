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
library(data.table) # Added for enhanced violin plot
library(shinycssloaders) # <--- ADD THIS LINE


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
        tabPanel("Feature Plot", plotlyOutput("featPlot", height = "600px")),
        tabPanel("Bubble Plot",  plotlyOutput("dotPlot",   height = "800px")),
        tabPanel("Violin Plot",  withSpinner(uiOutput("vlnPlot_ui"), type = 6)), # Updated for dynamic height and spinner
        tabPanel("Heatmap",      plotlyOutput("heatmapPlot",height = "800px")),
        tabPanel("Box Plot",     plotlyOutput("boxPlot",   height = "800px"))
      )
    )
  )
)

# --- SERVER LOGIC -----------------------------------------------------------
server <- function(input, output, session) {
  
  # 1) Load and merge expression + metadata -------------------------------
  # This reactive now primarily loads data for plots that don't need the full data.table processing
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
  
  # Data for the enhanced Violin Plot (using data.table for efficiency)
  expr_dt_full <- reactive({
    req(input$expr)
    dt <- fread(input$expr$datapath, encoding = "UTF-8")
    setnames(dt, names(dt)[1], "gene")
    # Clean gene names once when loading the file
    dt[, gene := toupper(trimws(gene))]
    dt
  })
  
  meta_dt_full <- reactive({
    req(input$meta)
    dt <- fread(input$meta$datapath)
    setnames(dt, names(dt)[1], "cell")
    dt[, cell := trimws(as.character(cell))]
    dt
  })
  
  # 2) Dynamic UI for gene, group, split ----------------------------------
  output$gene_ui <- renderUI({
    req(expr_dt_full()) # Use data.table reactive for gene choices
    gene_choices <- unique(expr_dt_full()$gene)
    selectizeInput("genes", "Select Gene(s):",
                   choices  = NULL, # Set to NULL and update via observeEvent
                   selected = NULL,
                   multiple = TRUE,
                   options = list(placeholder = "Type to search...", maxOptions = 100))
  })
  
  # Update gene choices using observeEvent for selectizeInput
  observeEvent(expr_dt_full(), {
    updateSelectizeInput(session, "genes",
                         choices = unique(expr_dt_full()$gene),
                         selected = NULL,
                         server = TRUE)
  })
  
  output$group_ui <- renderUI({
    req(meta_dt_full()) # Use data.table reactive for metadata columns
    meta_cols <- setdiff(names(meta_dt_full()), "cell")
    selectInput("group", "Group By (Categories/Colors):",
                choices  = meta_cols,
                selected = meta_cols[1])
  })
  output$split_ui <- renderUI({
    req(meta_dt_full()) # Use data.table reactive for metadata columns
    meta_cols <- setdiff(names(meta_dt_full()), "cell")
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
  
  # === VIOLIN PLOT SECTION START (Updated) ===
  # Dynamic plot container for Violin Plot
  output$vlnPlot_ui <- renderUI({
    req(input$genes)
    n    <- length(input$genes)
    # Adjust columns based on number of genes, max 2 columns for better readability
    cols <- if (n <= 2) n else 2
    rows <- ceiling(n/cols)
    height_px <- paste0(400*rows, "px") # Adjust height dynamically
    plotlyOutput("vlnPlot", height = height_px, width = "100%")
  })
  
  # Build and render violin plot
  output$vlnPlot <- renderPlotly({
    req(input$genes, input$group, expr_dt_full(), meta_dt_full())
    
    genes_selected <- toupper(trimws(input$genes)) # Clean selected gene names
    grp_col <- input$group
    spl_col <- input$split
    n_genes <- length(genes_selected)
    ncol_wrap <- if (n_genes <= 2) n_genes else 2
    
    withProgress(message = "Building violin plot(s)…", value = 0, {
      incProgress(0.1, detail = "Filtering selected genes")
      dt_s <- expr_dt_full()[gene %in% genes_selected]
      
      if (nrow(dt_s) == 0) {
        showNotification("No matching genes found. Check typos or formatting.", type = "error")
        return(NULL)
      }
      
      incProgress(0.2, detail = "Melting and joining data")
      df_long <- melt(dt_s,
                      id.vars = "gene",
                      variable.name = "cell",
                      value.name = "expression",
                      variable.factor = FALSE)
      df_long[, cell := trimws(as.character(cell))] # Clean cell names
      
      df <- merge(df_long, meta_dt_full(), by = "cell", all = FALSE)
      
      if (nrow(df) == 0) {
        showNotification("No matching cells between expression and metadata.", type = "error")
        return(NULL)
      }
      
      df[, group_col := factor(get(grp_col))]
      if (spl_col != "None") {
        df[, split_col := factor(get(spl_col))]
      } else {
        df[, split_col := NULL]
      }
      
      incProgress(0.15, detail = "Filtering small groups")
      if (spl_col != "None") {
        cnts <- df[, .N, by = .(gene, group_col, split_col)]
        keep <- cnts[N >= 2, .(gene, group_col, split_col)]
        df_big <- merge(df, keep, by = c("gene","group_col","split_col"))
      } else {
        cnts <- df[, .N, by = .(gene, group_col)]
        keep <- cnts[N >= 2, .(gene, group_col)]
        df_big <- merge(df, keep, by = c("gene","group_col"))
      }
      
      if (nrow(df_big) == 0) {
        showNotification("No groups have ≥2 cells after filtering. Try different settings.", type = "error")
        return(NULL)
      }
      
      incProgress(0.35, detail = "Rendering plot")
      if (spl_col != "None") {
        p <- ggplot(df_big,
                    aes(x = group_col, y = expression,
                        fill = split_col, group = interaction(group_col, split_col))) +
          geom_violin(trim = TRUE, color = "black", alpha = 0.6) +
          geom_jitter(aes(color = split_col),
                      position = position_jitterdodge(0.9, 0.15),
                      size = 0.6, alpha = 0.6) +
          scale_fill_viridis_d(name = spl_col) +
          scale_color_viridis_d(name = spl_col)
      } else {
        p <- ggplot(df_big,
                    aes(x = group_col, y = expression,
                        fill = group_col, group = group_col)) +
          geom_violin(trim = TRUE, color = "black", alpha = 0.6) +
          geom_jitter(aes(color = group_col),
                      position = position_jitter(0.15),
                      size = 0.6, alpha = 0.6) +
          scale_fill_viridis_d(guide = "none") +
          scale_color_viridis_d(guide = "none")
      }
      
      p <- p +
        facet_wrap(~gene, scales = "free_y", ncol = ncol_wrap, drop = FALSE) +
        labs(x = NULL, y = "Normalized Expression") +
        theme_minimal() +
        theme(
          panel.spacing   = unit(1, "cm"),
          axis.text.x     = element_text(angle = 45, hjust = 1),
          strip.text      = element_text(size = 12),
          plot.margin     = margin(10, 10, 10, 10)
        )
      
      ggplotly(p, tooltip = c("group_col", "split_col", "expression")) %>%
        layout(
          margin = list(b = 150),
          legend = list(orientation = "h", y = -0.2, x = 0.5, xanchor = "center")
        )
    })
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
