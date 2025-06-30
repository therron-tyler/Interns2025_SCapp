library(shiny)
library(ggplot2)

# Empty placeholder plot function
emptyPlot <- function(title) {
  plot(1, type = "n", xlab = "", ylab = "", main = title, axes = FALSE)
  box()
}
# Sets max file size to 9 GB
options(shiny.maxRequestSize = 9437184000) 

# Tab names
plot_names <- c("Violin Plot", "Heatmap", "Feature Plot", "Box & Whiskers Plot", "Bubble Plot")

# UI
ui <- fluidPage(
  titlePanel("Winter Lab Single Cell Analysis"),
  
  # Two file inputs
  fileInput("expr_file", "Upload Expression Data (.tsv)", accept = ".tsv"),
  fileInput("meta_file", "Upload Metadata (.tsv)", accept = ".tsv"),
  
  # Tab panels (unchanged)
  do.call(tabsetPanel, lapply(seq_along(plot_names), function(i) {
    tabPanel(plot_names[i], plotOutput(paste0("plot", i)))
  }))
)

# Server
server <- function(input, output) {
  # Reactive: load expression data
  expr_data <- reactive({
    req(input$expr_file)
    read.delim(input$expr_file$datapath, header = TRUE, row.names = 1, check.names = FALSE)
  })
  
  # Reactive: load metadata
  meta_data <- reactive({
    req(input$meta_file)
    read.delim(input$meta_file$datapath, header = TRUE, row.names = 1, check.names = FALSE)
  })
  
  # Generate placeholder plots 
  for (i in seq_along(plot_names)) {
    local({
      index <- i
      output[[paste0("plot", index)]] <- renderPlot({
        emptyPlot(paste0(plot_names[index], " Placeholder"))
      })
    })
  }
}

shinyApp(ui = ui, server = server)