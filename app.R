#here, set shiny upload limit to 600 MB
options(shiny.maxRequestSize = 600 * 1024^2)

#this loads shiny package for building interactive web apps
library(shiny)
# loads Seurat package for single-cell analysis
library(Seurat)
#loads readr for reading TSV files
library(readr)
#loads tibble for tidy data frames
library(tibble)
#loads ggplot2 for plotting
library(ggplot2)
#loads pheatmap for heatmap visualization
library(pheatmap)

#this defines the user interface of the app
ui <- fluidPage(
  #this sets the application title
  titlePanel("Single-cell Explorer (Seurat + Shiny)"),
  #here, define the layout with sidebar and main area
  sidebarLayout(
    sidebarPanel(
      #here, input for expression matrix file
      fileInput("expr", "Expression matrix (.tsv or .tsv.gz)",
                accept = c(".tsv", ".tsv.gz")),
      #here, input for metadata file
      fileInput("meta", "Metadata (.tsv or .tsv.gz)",
                accept = c(".tsv", ".tsv.gz")),
      #this adds a horizontal rule for separation
      hr(),
      #here, dynamic UI placeholder for gene selection
      uiOutput("gene_ui"),
      #here, dynamic UI placeholder for grouping selection
      uiOutput("group_ui")
    ),
    mainPanel(
      #this creates tabs for different plot types
      tabsetPanel(
        tabPanel("Violin", plotOutput("vlnPlot")),
        tabPanel("Feature", plotOutput("featPlot")),
        tabPanel("Bubble", plotOutput("dotPlot")),
        tabPanel("Heatmap", plotOutput("heatmapPlot")),
        tabPanel("Boxplot", plotOutput("boxPlot"))
      )
    )
  )
)

#here, define server logic
server <- function(input, output, session) {
  #this reactive expression reads data and builds Seurat object
  sc <- reactive({
    #here, ensure expression and metadata inputs are available
    req(input$expr, input$meta)
    #this reads the expression matrix
    expr_df <- read_tsv(input$expr$datapath)
    #this extracts gene names from first column
    genes <- expr_df[[1]]
    #this removes the first column to keep numeric data
    expr_df <- expr_df[,-1]
    #this converts data frame to matrix
    mat <- as.matrix(expr_df)
    #here, assign gene names as row names
    rownames(mat) <- genes
    
    #this reads metadata file
    meta_df <- read_tsv(input$meta$datapath)
    #this assumes first column is cell barcode
    rownames(meta_df) <- meta_df[[1]]
    #here, drop the barcode column after setting row names
    meta_df <- meta_df[,-1]
    
    #this creates a Seurat object from counts and metadata
    seu <- CreateSeuratObject(counts = mat, meta.data = meta_df)
    
    #now i'll re-read metadata to preserve row names directly
    meta_df <- read.delim(
      input$meta$datapath,
      sep = "\t",
      header = TRUE,
      row.names = 1,
      stringsAsFactors = FALSE
    )
    
    #this extracts real UMAP coordinates from metadata
    cells <- Cells(seu)
    umap_coords <- meta_df[cells, c("UMAP_Xaxis","UMAP_Yaxis"), drop = FALSE]
    
    #here, rename columns to match Seurat format
    colnames(umap_coords) <- paste0("UMAP_", seq_len(ncol(umap_coords)))
    
    #this creates DimReduc object for UMAP embedding
    seu[["umap"]] <- CreateDimReducObject(
      embeddings = as.matrix(umap_coords),
      key = "UMAP_",
      assay = DefaultAssay(seu)
    )
    
    #this normalizes data
    seu <- NormalizeData(seu, verbose = FALSE)
    #this identifies variable features
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    #this scales data across all genes
    seu <- ScaleData(seu, features = rownames(seu), verbose = FALSE)
    
    #this returns the Seurat object with all preprocessing done
    seu
  })
  
  #this dynamically populates gene selection based on sc() object
  output$gene_ui <- renderUI({
    genes <- rownames(sc())
    selectInput("genes", "Select gene(s):", choices = genes,
                selected = genes[1], multiple = TRUE)
  })
  #this dynamically populates grouping options based on metadata columns
  output$group_ui <- renderUI({
    md <- sc()[[]]
    groups <- colnames(md)
    selectInput("group", "Group by metadata:", choices = groups,
                selected = "seurat_clusters")
  })
  
  #this renders violin plot of gene expression by group
  output$vlnPlot <- renderPlot({
    req(input$genes, input$group)
    VlnPlot(sc(), features = input$genes,
            group.by = input$group, pt.size = 0.1)
  })
  
  #this renders feature plot using UMAP reduction
  output$featPlot <- renderPlot({
    req(input$genes)
    FeaturePlot(sc(), features = input$genes, reduction = "umap")
  })
  
  #this renders bubble/dot plot of expression across groups
  output$dotPlot <- renderPlot({
    req(input$genes, input$group)
    DotPlot(sc(), features = input$genes,
            group.by = input$group) + RotatedAxis()
  })
  
  #this renders heatmap of selected genes
  output$heatmapPlot <- renderPlot({
    req(input$genes)
    #this checks if at least two genes are selected
    if (length(input$genes) < 2) {
      plot.new()
      text(0.5, 0.5, "Select at least two genes for a heatmap", cex = 1.2)
      return()
    }
    #this fetches scaled data for selected genes
    mat <- GetAssayData(
      object = sc(),
      assay = "RNA",
      slot = "scale.data"
    )[input$genes, , drop = FALSE]
    #this ensures data is a plain matrix
    mat <- as.matrix(mat)
    #here, plot heatmap without additional scaling
    pheatmap(
      mat,
      scale = "none",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = TRUE,
      show_colnames = FALSE
    )
  })
  
  #this renders boxplot of gene expression by group
  output$boxPlot <- renderPlot({
    req(input$genes, input$group)
    df <- FetchData(sc(), vars = c(input$genes, input$group))
    colnames(df) <- c("expression", "group")
    #this creates ggplot boxplot with rotated axis labels
    ggplot(df, aes(x = group, y = expression)) +
      geom_boxplot() +
      xlab(input$group) + ylab("Expression") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
}

#this launches the Shiny app with defined UI and server
shinyApp(ui, server)
