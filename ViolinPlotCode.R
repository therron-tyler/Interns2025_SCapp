#–––– Increase upload limit –––––
options(shiny.maxRequestSize = 600 * 1024^2)

#–––– Packages –––––
# install.packages(c("shiny","data.table","ggplot2","viridis","shinycssloaders"))
library(shiny)
library(data.table)      # fread, melt
library(ggplot2)         # ggplot2 + position_jitter / position_jitterdodge
library(viridis)         # scale_fill_viridis_d()
library(shinycssloaders) # withSpinner()

#–––– UI –––––
ui <- fluidPage(
  titlePanel("Single-cell Violin Plot Explorer"),
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
      withSpinner(plotOutput("vlnPlot", width="100%"))
    )
  )
)

#–––– Server –––––
server <- function(input, output, session) {
  
  # 1) Read header row for real column names
  expr_header <- reactive({
    req(input$expr)
    names(fread(input$expr$datapath, nrows = 0))
  })
  
  # 2) Read only first column (genes) to populate dropdown
  expr_genes <- reactive({
    req(input$expr)
    fread(input$expr$datapath, select = 1)[[1]]
  })
  
  # 3) Load metadata once
  meta_dt <- reactive({
    req(input$meta)
    dt <- fread(input$meta$datapath)
    setnames(dt, names(dt)[1], "cell")
    dt
  })
  
  # 4) Populate selectors immediately
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
  
  # 5) Build & draw the violin plot
  output$vlnPlot <- renderPlot({
    req(input$expr, input$meta, input$genes, input$group)
    
    # fundamentals
    path      <- input$expr$datapath
    genes     <- input$genes
    hdr       <- expr_header()
    meta      <- meta_dt()
    grp       <- input$group
    spl       <- input$split_by
    n_genes   <- length(genes)
    ncol_wrap <- if (n_genes <= 2) n_genes else 2
    
    withProgress(message = "Building violin(s)…", value = 0, {
      # 5.1 build grep pattern
      incProgress(0.1, detail = "Building grep pattern")
      esc     <- gsub("([\\^\\$\\.\\|\\(\\)\\[\\]\\*\\+\\?\\\\])","\\\\\\\\1", genes)
      pattern <- paste0("^(", paste(esc, collapse="|"), ")\t")
      
      # 5.2 extract matching lines
      incProgress(0.1, detail = "Extracting lines")
      cmd      <- paste("grep -E", shQuote(pattern), shQuote(path))
      dt_s     <- tryCatch(fread(cmd = cmd, header = FALSE, sep = "\t"), error = function(e) NULL)
      if (is.null(dt_s) || nrow(dt_s)==0) {
        showNotification("No matching genes found.", type="error")
        return(NULL)
      }
      
      # 5.3 reassign column names
      incProgress(0.1, detail = "Reassigning column names")
      setnames(dt_s, hdr)            # real header
      setnames(dt_s, hdr[1], "gene") # first col → "gene"
      
      # 5.4 melt to long
      incProgress(0.15, detail = "Melting to long form")
      df_long <- melt(
        dt_s,
        id.vars         = "gene",
        variable.name   = "cell",
        value.name      = "expression",
        variable.factor = FALSE
      )
      
      # 5.5 join metadata
      incProgress(0.1, detail = "Joining metadata")
      df <- meta[df_long, on = "cell", nomatch = 0]
      df[, group := factor(get(grp))]
      if (nzchar(spl)) df[, split := factor(get(spl))] else df[, split := NULL]
      
      # 5.6 prefilter for valid violin density
      incProgress(0.15, detail = "Filtering groups for violin")
      if (nzchar(spl)) {
        cnts   <- df[, .N, by = .(gene, group, split)]
        keep   <- cnts[N >= 2, .(gene, group, split)]
        df_big <- merge(df, keep, by = c("gene","group","split"))
      } else {
        cnts   <- df[, .N, by = .(gene, group)]
        keep   <- cnts[N >= 2, .(gene, group)]
        df_big <- merge(df, keep, by = c("gene","group"))
      }
      
      # 5.7 render
      incProgress(0.35, detail = "Rendering plot")
      if (nzchar(spl)) {
        p <- ggplot() +
          geom_violin(data = df_big,
                      aes(x = group, y = expression, fill = split,
                          group = interaction(group,split)),
                      trim = TRUE, color = "black", alpha = 0.6) +
          geom_jitter(data = df,
                      aes(x = group, y = expression, color = split,
                          group = interaction(group,split)),
                      position = position_jitterdodge(dodge.width = 0.9,
                                                      jitter.width = 0.15),
                      size = 0.6, alpha = 0.6) +
          scale_fill_viridis_d(name = spl) +
          scale_color_viridis_d(name = spl)
      } else {
        p <- ggplot() +
          geom_violin(data = df_big,
                      aes(x = group, y = expression, fill = group, group = group),
                      trim = TRUE, color = "black", alpha = 0.6) +
          geom_jitter(data = df,
                      aes(x = group, y = expression, color = group, group = group),
                      position = position_jitter(width = 0.15),
                      size = 0.6, alpha = 0.6) +
          scale_fill_viridis_d(guide = "none") +
          scale_color_viridis_d(guide = "none")
      }
      
      print(
        p +
          facet_wrap(~ gene, scales = "free_y", ncol = ncol_wrap, drop = FALSE) +
          labs(x = NULL, y = "Normalized Expression") +
          theme_minimal() +
          theme(
            panel.spacing = unit(1, "cm"),
            axis.text.x   = element_text(angle = 45, hjust = 1),
            strip.text    = element_text(size = 12),
            plot.margin   = margin(10, 10, 10, 10)
          )
      )
    })
  },
  # dynamic height: 350px per row
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

