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
library(patchwork) # Used for combining plots
library(RColorBrewer)
library(data.table)
library(shinycssloaders)
library(htmlwidgets) # <-- ADD THIS LINE for saving HTML plots
library(colorspace)


# Increase upload limit
options(shiny.maxRequestSize = 600 * 1024^2)
Sys.setenv("VROOM_CONNECTION_SIZE" = 10000000)

# --- UI DEFINITION -----------------------------------------------------------
ui <- fluidPage(
  titlePanel("Interactive Gene Expression Visualization Tool - Sohum, Suki, Kim, and Zach"),
  sidebarLayout(
    sidebarPanel(
      fileInput("expr", "1. Expression Matrix (.tsv or .tsv.gz)", accept = c(".tsv", ".tsv.gz")),
      fileInput("meta", "2. Metadata File (.tsv or .tsv.gz)",   accept = c(".tsv", ".tsv.gz")),
      hr(),
      uiOutput("gene_ui"),      # Gene selector
      uiOutput("group_ui"),     # Group-by selector
      uiOutput("split_ui"),     # Split-by selector
      checkboxInput("color_by_gene",
                    "Feature Plot: Color by gene expression",
                    value = FALSE),
      p(strong("Note:"), "For Violin/Box plots, multiple genes will be shown in a grid."),
      
      # --- NEW: Download Center UI ---
      hr(),
      wellPanel(
        h4("Download Center"),
        p("Generate and download the combined UMAP feature plots below."),
        downloadButton("downloadPng", "Download as PNG", icon = icon("download")),
        downloadButton("downloadHtml", "Download as HTML", icon = icon("download"))
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Feature Plot", withSpinner(plotlyOutput("featPlot", height = "900px"), type = 6)),
        tabPanel("Bubble Plot",  plotlyOutput("dotPlot",     height = "800px")),
        tabPanel("Violin Plot",  withSpinner(uiOutput("vlnPlot_ui"), type = 6)),
        tabPanel("Heatmap",      plotlyOutput("heatmapPlot", height = "800px")),
        tabPanel("Box Plot",     plotlyOutput("boxPlot",     height = "800px"))
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

  expr_dt_full <- reactive({
    req(input$expr)
    dt <- fread(input$expr$datapath, encoding = "UTF-8")
    setnames(dt, names(dt)[1], "gene")
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
    req(expr_dt_full())
    gene_choices <- unique(expr_dt_full()$gene)
    selectizeInput("genes", "Select Gene(s):",
                   choices  = NULL,
                   selected = NULL,
                   multiple = TRUE,
                   options = list(placeholder = "Type to search...", maxOptions = 100))
  })

  observeEvent(expr_dt_full(), {
    updateSelectizeInput(session, "genes",
                         choices = unique(expr_dt_full()$gene),
                         selected = NULL,
                         server = TRUE)
  })

  output$group_ui <- renderUI({
    req(meta_dt_full())
    meta_cols <- setdiff(names(meta_dt_full()), "cell")
    selectInput("group", "Group By (Categories/Colors):",
                choices  = meta_cols,
                selected = meta_cols[1])
  })
  output$split_ui <- renderUI({
    req(meta_dt_full())
    meta_cols <- setdiff(names(meta_dt_full()), "cell")
    selectInput("split", "Split By (Rows/Panels):",
                choices  = c("None", meta_cols),
                selected = "None")
  })

  # === FEATURE PLOT SECTION START (RESTRUCTURED FOR DOWNLOADS) ===
  
  # --- NEW: Reactive ggplot object for sharing ---
  # This reactive builds the ggplot object, which can then be used by both
  # the renderPlotly output and the download handlers.
  feature_plot_object <- reactive({
      # Ensure base data and UMAP coordinates are available
      req(plot_data(), "UMAP_Xaxis" %in% names(plot_data()), input$group)
      df <- plot_data()

      # Prepare enhanced hover text based on user selections
      if (!is.null(input$genes) && length(input$genes) > 0) {
          if (length(input$genes) == 1) { # Case for a SINGLE selected gene
              gene <- input$genes[1]
              summary_stats <- df %>%
                  group_by(.data[[input$group]]) %>%
                  summarise(
                      avg_expr = mean(.data[[gene]], na.rm = TRUE),
                      pct_expr = 100 * sum(.data[[gene]] > 0, na.rm = TRUE) / n(),
                      .groups = "drop"
                  )
              df <- left_join(df, summary_stats, by = input$group)
              df$hover_text <- paste0(
                  "<b>Cell:</b> ", df$cell, "<br>",
                  "<b>Group:</b> ", df[[input$group]], "<br>",
                  #"<hr>", # Adds a horizontal line
                  "<b>Gene:</b> ", gene, "<br>",
                  "<b>Expression (this cell):</b> ", round(df[[gene]], 2), "<br>",
                  "<b>Avg. Expr (group):</b> ", round(df$avg_expr, 2), "<br>",
                  "<b>% Expressing (group):</b> ", round(df$pct_expr, 1), "%"
              )
          } else { # Case for MULTIPLE selected genes (Module Score)
              df$ModuleScore <- rowMeans(df[, input$genes, drop = FALSE], na.rm = TRUE)
              summary_stats <- df %>%
                  group_by(.data[[input$group]]) %>%
                  summarise(avg_score = mean(ModuleScore, na.rm = TRUE), .groups = "drop")
              df <- left_join(df, summary_stats, by = input$group)
              df$hover_text <- paste0(
                  "<b>Cell:</b> ", df$cell, "<br>",
                  "<b>Group:</b> ", df[[input$group]], "<br>",
                  "<hr>",
                  "<b>Module Score (this cell):</b> ", round(df$ModuleScore, 2), "<br>",
                  "<b>Avg. Score (group):</b> ", round(df$avg_score, 2)
              )
          }
      } else { # Case for NO genes selected
          df$hover_text <- paste0("<b>Cell:</b> ", df$cell, "<br><b>Group:</b> ", df[[input$group]])
      }


      # --- Case 1: Color by Metadata ---
      if (!input$color_by_gene) {
          p <- ggplot(df, aes(x = UMAP_Xaxis, y = UMAP_Yaxis, color = .data[[input$group]], text = hover_text)) +
              geom_point(size = 1, alpha = 0.8) +
              guides(color = guide_legend(override.aes = list(size = 3))) +
              labs(title = "UMAP Colored by Metadata Group", color = input$group) +
              theme_minimal() +
              theme(legend.title = element_text(size = 10), legend.text = element_text(size = 8))
          return(p)
      }

      # --- Case 2: Color by Gene Expression ---
      req(input$genes)

      # --- Subcase 2a: NO Splitting ---
      if (input$split == "None") {
          if (length(input$genes) == 1) {
              p <- ggplot(df, aes(x = UMAP_Xaxis, y = UMAP_Yaxis, color = .data[[input$genes]], text = hover_text)) +
                  geom_point(size = 1, alpha = 0.8) +
                  scale_color_gradient(low = "lightgrey", high = "darkviolet") +
                  labs(title = paste("Expression of", input$genes), color = "Expression") +
                  theme_minimal()
          } else {
              p <- ggplot(df, aes(x = UMAP_Xaxis, y = UMAP_Yaxis, color = ModuleScore, text = hover_text)) +
                  geom_point(size = 1, alpha = 0.8) +
                  scale_color_gradient(low = "lightgrey", high = "darkviolet") +
                  labs(title = "Module Score for Selected Genes", color = "Avg. Expr.") +
                  theme_minimal()
          }
          return(p)
      }

      # --- Subcase 2b: WITH Splitting ---
      p_meta <- ggplot(df, aes(x = UMAP_Xaxis, y = UMAP_Yaxis, color = .data[[input$group]], text = hover_text)) +
          geom_point(size = 1.5, alpha = 0.9) +
          guides(color = guide_legend(override.aes = list(size = 4))) +
          labs(title = paste("UMAP grouped by", input$group), color = input$group) +
          theme_minimal()

      df_long <- df %>%
          select(UMAP_Xaxis, UMAP_Yaxis, cell, all_of(input$group), all_of(input$split), all_of(input$genes)) %>%
          pivot_longer(
              cols = all_of(input$genes),
              names_to = "gene",
              values_to = "expression"
          )

      summary_split <- df_long %>%
          group_by(gene, .data[[input$split]]) %>%
          summarise(
              avg_expr_split = mean(expression, na.rm = TRUE),
              pct_expr_split = 100 * sum(expression > 0, na.rm = TRUE) / n(),
              .groups = "drop"
          )
      df_long <- left_join(df_long, summary_split, by = c("gene", input$split))
      
      df_long$hover_text_split <- paste0(
          "<b>Cell:</b> ", df_long$cell, "<br>",
          "<b>Group:</b> ", df_long[[input$group]], "<br>",
          "<b>", input$split, ":</b> ", df_long[[input$split]], "<br>",
          "<hr>",
          "<b>Gene:</b> ", df_long$gene, "<br>",
          "<b>Expression (this cell):</b> ", round(df_long$expression, 2), "<br>",
          "<b>Avg. Expr (this panel):</b> ", round(df_long$avg_expr_split, 2), "<br>",
          "<b>% Expressing (this panel):</b> ", round(df_long$pct_expr_split, 1), "%"
      )

      df_bg <- df_long %>% distinct(UMAP_Xaxis, UMAP_Yaxis, .data[[input$split]], gene)

      # --- FIX: Moved the `text` aesthetic from the main ggplot call to the specific geom_point layer ---
      p_split_feat <- ggplot(df_long, aes(x = UMAP_Xaxis, y = UMAP_Yaxis)) +
          geom_point(data = df_bg, color = "grey90", size = 0.5, alpha = 0.5) +
          geom_point(data = filter(df_long, expression > 0), aes(color = expression, text = hover_text_split), size = 0.8) +
          scale_color_gradient(low = "lightgrey", high = "darkviolet") +
          facet_grid(
              rows = vars(.data[[input$split]]),
              cols = vars(gene)
          ) +
          labs(x = NULL, y = NULL, color = "Expression") +
          theme_void() +
          theme(
              strip.text.x = element_text(size = 12, face = "bold"),
              strip.text.y = element_text(size = 12, face = "bold", angle = 0),
              panel.border = element_rect(color = "black", fill = NA, size = 0.5),
              legend.position = "bottom"
          )

      combined_plot <- p_meta / p_split_feat +
                       plot_layout(heights = c(1, length(unique(df[[input$split]]))))

      return(combined_plot)
  })

  # Render the plotly object for the UI
  output$featPlot <- renderPlotly({
      plot_gg <- feature_plot_object()
      req(plot_gg)
      ggplotly(plot_gg, tooltip = "text")
  })
  
  # --- NEW: Download Handlers ---
  output$downloadPng <- downloadHandler(
    filename = function() {
      paste0("feature_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      plot_to_save <- feature_plot_object()
      req(plot_to_save)
      # Use ggsave for high-quality export. Adjust dimensions as needed.
      ggsave(file, plot = plot_to_save, device = "png", width = 10, height = 12, dpi = 300, units = "in")
    }
  )
  
  output$downloadHtml <- downloadHandler(
    filename = function() {
      paste0("feature_plot_", Sys.Date(), ".html")
    },
    content = function(file) {
      plot_gg <- feature_plot_object()
      req(plot_gg)
      # Convert to plotly and then save as a self-contained HTML file
      plot_ly <- ggplotly(plot_gg, tooltip = "text")
      saveWidget(widget = plot_ly, file = file, selfcontained = TRUE)
    }
  )
  
  # === FEATURE PLOT SECTION END ===

  # === BUBBLE PLOT SECTION START ===
  output$dotPlot <- renderPlotly({
    req(input$genes, input$group, plot_data())
    df <- plot_data()
    
    split_var <- if (input$split != "None") input$split else NULL
    group_var <- input$group
    
    dot_data <- df %>%
      select(all_of(c(input$genes, split_var, group_var))) %>%
      pivot_longer(cols = all_of(input$genes), names_to = "gene", values_to = "expression")
    
    summary_df <- dot_data %>%
      group_by(gene, .data[[split_var]], .add = TRUE) %>%
      { if (!is.null(group_var)) group_by(., .data[[group_var]], .add = TRUE) else . } %>%
      summarise(
        avg_expr = mean(expression, na.rm = TRUE),
        pct_expr = 100 * sum(expression > 0, na.rm = TRUE) / n(),
        .groups = "drop"
      )
    
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
    
    if (!is.null(group_var)) {
      p <- p + facet_wrap(vars(.data[[group_var]]), ncol = 1)
    }
    
    ggplotly(p, tooltip = "text") %>% layout(legend = list(orientation = "h", y = -0.2))
  })
  # === BUBBLE PLOT SECTION END ===

  # === VIOLIN PLOT SECTION START (Enhanced) ===
  
  # Helper function for better colors
  get_colors <- function(n) {
    if (n <= 8) {
      brewer.pal(n, "Pastel1")
    } else if (n <= 12) {
      brewer.pal(n, "Set3")
    } else {
      qualitative_hcl(n, palette = "Pastel 1")
    }
  }
  
  output$vlnPlot_ui <- renderUI({
    req(input$genes)
    n    <- length(input$genes)
    cols <- if (n <= 2) n else 2
    rows <- ceiling(n/cols)
    height_px <- paste0(400*rows, "px")
    plotlyOutput("vlnPlot", height = height_px, width = "100%")
  })

  output$vlnPlot <- renderPlotly({
    req(input$genes, input$group, expr_dt_full(), meta_dt_full())
    
    genes_selected <- toupper(trimws(input$genes))
    grp_col <- input$group
    spl_col <- input$split
    n_genes <- length(genes_selected)
    ncol_wrap <- if (n_genes <= 2) n_genes else 2
    
    withProgress(message = "Building violin plot(s)…", value = 0, {
      incProgress(0.25, detail = "Reading expression data...")
      dt_s <- expr_dt_full()[gene %in% genes_selected]
      
      if (nrow(dt_s) == 0) {
        showNotification("Selected genes not found in the expression matrix.", type = "error")
        return(NULL)
      }
      
      incProgress(0.1, detail = "Melting to long form")
      df_long <- melt(dt_s,
                      id.vars = "gene",
                      variable.name = "cell",
                      value.name = "expression",
                      variable.factor = FALSE)
      df_long[, cell := trimws(as.character(cell))]
      
      incProgress(0.1, detail = "Joining metadata")
      df <- merge(df_long, meta_dt_full(), by = "cell", all = FALSE)
      
      if (nrow(df) == 0) {
        showNotification("No matching cells between expression and metadata.", type = "error")
        return(NULL)
      }
      
      df[, group := factor(get(grp_col))]
      if (spl_col != "None" && nzchar(spl_col)) {
        df[, split := factor(get(spl_col))]
      } else {
        df[, split := NULL]
      }
      
      # Enhanced tooltip creation with metadata
      meta_cols <- setdiff(names(meta_dt_full()), "cell")
      df[, tooltip_text := {
        meta_info <- sapply(meta_cols, function(col) paste0(col, ": ", get(col)))
        paste(c(
          paste("Cell:", cell),
          paste0("Expression: ", round(expression, 3)),
          meta_info
        ), collapse = "<br>")
      }, by = seq_len(nrow(df))]
      
      incProgress(0.1, detail = "Filtering for densities")
      if (spl_col != "None" && nzchar(spl_col)) {
        cnts <- df[, .N, by = .(gene, group, split)]
        keep <- cnts[N >= 2, .(gene, group, split)]
        df_big <- merge(df, keep, by = c("gene", "group", "split"))
      } else {
        cnts <- df[, .N, by = .(gene, group)]
        keep <- cnts[N >= 2, .(gene, group)]
        df_big <- merge(df, keep, by = c("gene", "group"))
      }
      
      if (nrow(df_big) == 0) {
        showNotification("No groups have ≥2 cells after filtering. Try different settings.", type = "error")
        return(NULL)
      }
      
      incProgress(0.2, detail = "Rendering plot")
      
      # Generate appropriate number of distinct pastel colors
      n_groups <- if (spl_col != "None" && nzchar(spl_col)) {
        length(unique(df_big$split))
      } else {
        length(unique(df_big$group))
      }
      pastel_colors <- get_colors(n_groups)
      
      if (spl_col != "None" && nzchar(spl_col)) {
        p <- ggplot(df_big,
                    aes(x = group, y = expression,
                        fill = split, group = interaction(group, split),
                        text = tooltip_text)) +
          geom_violin(trim = TRUE, color = "black", alpha = 0.6) +
          geom_jitter(aes(color = split),
                      position = position_jitterdodge(0.9, 0.15),
                      size = 0.6, alpha = 0.6) +
          scale_fill_manual(name = spl_col, values = pastel_colors) +
          scale_color_manual(values = pastel_colors, guide = "none")
      } else {
        p <- ggplot(df_big,
                    aes(x = group, y = expression,
                        fill = group, group = group,
                        text = tooltip_text)) +
          geom_violin(trim = TRUE, color = "black", alpha = 0.6) +
          geom_jitter(aes(color = group),
                      position = position_jitter(0.15),
                      size = 0.6, alpha = 0.6) +
          scale_fill_manual(values = pastel_colors, guide = "none") +
          scale_color_manual(values = pastel_colors, guide = "none")
      }
      
      p <- p +
        facet_wrap(~ gene, scales = "free_y",
                   ncol = ncol_wrap, drop = FALSE) +
        labs(x = NULL, y = "Normalized Expression") +
        theme_minimal() +
        theme(
          panel.spacing = unit(1, "cm"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(size = 12),
          plot.margin = margin(10, 10, 10, 10)
        )
      
      ggplotly(p, tooltip = "text") %>%
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
    
    if (input$split != "None") {
      grp <- paste(df[[input$group]], df[[input$split]], sep = "_")
    } else {
      grp <- df[[input$group]]
    }
    
    mat <- as.matrix(df[input$genes])
    avg_mat <- t(rowsum(mat, grp) / as.vector(table(grp)))
    
    mat_z <- t(scale(t(avg_mat)))
    mat_z[is.na(mat_z)] <- 0
    
    plot_ly(
      x = colnames(mat_z),
      y = rownames(mat_z),
      z = mat_z,
      type = "heatmap",
      colors = colorRamp(c("#440154", "white", "#21918c")),
      showscale = TRUE
    ) %>%
      layout(
        title = list(
          text = "Expression Heatmap<br><sub>Each group = CellType_DiseaseStatus</sub>",
          x = 0.5,
          xanchor = "center"
        ),
        xaxis = list(title = "Group", tickangle = 45),
        yaxis = list(title = "Gene"),
        margin = list(l = 80, r = 20, b = 100, t = 80)
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
