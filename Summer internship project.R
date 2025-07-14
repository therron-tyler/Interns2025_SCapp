#–––– Increase upload limit –––––
options(shiny.maxRequestSize = 600 * 1024^2)

#–––– Packages –––––
library(shiny)
library(data.table)
library(ggplot2)
library(viridis)
library(shinycssloaders)

#–––– UI –––––
ui <- fluidPage(
  titlePanel("Single-cell Feature Plot Explorer"),  # Changed title
  sidebarLayout(
    sidebarPanel(
      fileInput("expr",   "Expression matrix (.tsv)", accept = ".tsv"),
      fileInput("meta",   "Metadata (.tsv)",        accept = ".tsv"),
      hr(),
      selectizeInput(
        "genes", "Select gene(s):",
        choices  = NULL, multiple = TRUE,
        options  = list(placeholder="Type to search…", maxOptions=100)
      ),
      selectInput("group",    "Group by (required):",    choices = NULL),
      selectInput("split_by", "Split by (optional):",    choices = NULL)
    ),
    mainPanel(
      withSpinner(plotOutput("vlnPlot", width="100%"))  # Kept ID as "vlnPlot" for minimal changes
    )
  )
)

#–––– Server –––––
server <- function(input, output, session) {
  
  expr_header <- reactive({
    req(input$expr)
    names(fread(input$expr$datapath, nrows = 0))
  })
  
  expr_genes <- reactive({
    req(input$expr)
    fread(input$expr$datapath, select = 1)[[1]]
  })
  
  meta_dt <- reactive({
    req(input$meta)
    dt <- fread(input$meta$datapath)
    setnames(dt, names(dt)[1], "cell")
    dt
  })
  
  observeEvent(expr_genes(), {
    genes <- expr_genes()
    updateSelectizeInput(session, "genes",
                         choices  = genes,
                         selected = genes[1],
                         server   = TRUE)
  })
  observeEvent(meta_dt(), {
    cols <- setdiff(names(meta_dt()), "cell")
    updateSelectInput(session, "group",    choices = cols,      selected = cols[1])
    updateSelectInput(session, "split_by", choices = c(None = "", cols), selected = "")
  })
  
  output$vlnPlot <- renderPlot({
    req(input$expr, input$meta, input$genes)
    path      <- input$expr$datapath
    genes     <- input$genes
    hdr       <- expr_header()
    meta      <- meta_dt()
    n_genes   <- length(genes)
    ncol_wrap <- if (n_genes <= 2) n_genes else 2
    
    withProgress(message = "Building Feature Plot(s)…", value = 0, {
      esc     <- gsub("([\\^\\$\\.\\|\\(\\)\\[\\]\\*\\+\\?\\\\])","\\\\\\\\1", genes)
      pattern <- paste0("^(", paste(esc, collapse="|"), ")\t")
      cmd     <- paste("grep -E", shQuote(pattern), shQuote(path))
      dt_s    <- tryCatch(fread(cmd = cmd, header = FALSE, sep = "\t"), error = function(e) NULL)
      if (is.null(dt_s) || nrow(dt_s)==0) {
        showNotification("No matching genes found.", type="error")
        return(NULL)
      }
      
      setnames(dt_s, hdr)
      setnames(dt_s, hdr[1], "gene")
      
      df_long <- melt(dt_s,
                      id.vars = "gene",
                      variable.name = "cell",
                      value.name = "expression",
                      variable.factor = FALSE)
      
      df <- meta[df_long, on = "cell", nomatch = 0]
      
      if (!all(c("UMAP_X", "UMAP_Y") %in% names(df))) {
        showNotification("Metadata must contain UMAP_X and UMAP_Y columns.", type = "error")
        return(NULL)
      }
      
      ggplot(df, aes(x = UMAP_X, y = UMAP_Y, color = expression)) +
        geom_point(alpha = 0.7, size = 1.2) +
        scale_color_viridis(option = "D") +
        facet_wrap(~ gene, scales = "free", ncol = ncol_wrap) +
        theme_minimal() +
        labs(x = "UMAP X", y = "UMAP Y", color = "Expression") +
        theme(
          strip.text = element_text(size = 12),
          axis.text = element_text(size = 10),
          plot.margin = margin(10, 10, 10, 10)
        )
    })
  },
  height = function() {
    req(input$genes)
    n    <- length(input$genes)
    cols <- if (n <= 2) n else 2
    rows <- ceiling(n / cols)
    350 * rows
  })
}

#–––– Launch App –––––
shinyApp(ui, server)
