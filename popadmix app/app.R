library(shiny)
source("appfunction.R", local = TRUE)
require(popadmix)
require(RColorBrewer)
require(ggplot2)
require(ggtern)

ui <- fluidPage(
  sidebarPanel(textInput("allele", "Enterloci"),
               checkboxGroupInput("choose", "Please choose", choices = c("Hadm", "Fst"), selected = NULL),
               
               
               selectInput("Breaks", "Set Breaks", c("Default", "Upload Break file")),
               uiOutput("ui1"),
               selectInput("Labels", "Set Labels", c("Default", "Upload Label file")),
               uiOutput("ui2"),
               
               tags$p("Break values:"),
               verbatimTextOutput("dynamic_value1"),
               
               tags$p("Label values:"),
               verbatimTextOutput("dynamic_value2"),
               
               actionButton("BUTTON1","Press the button")),
  
  mainPanel(plotOutput("Plot", width = "100%", height = "400px"))
)


server <- function(input, output, session){
  
  output$ui1 <- renderUI({
    switch(input$Breaks, "Upload Break file" = fileInput("file1", "Upload Break File"))
  })
  
  output$ui2 <- renderUI({
    switch(input$Labels, "Upload Label file" = fileInput("file2", "Upload Label File"))
  })
  
  
  observeEvent(
    input$BUTTON1, 
    {
      output$Plot <- renderPlot({
        if(input$choose == "Hadm"){
          if(is.null(input$file1)){new.break.file = NULL
          if(is.null(input$file2)){new.labels.file = NULL}
          plot_HorF(input$allele, y = input$choose)}  # file1 should be a .txt file
          else{
            ipfile1 <- input$file1
            break.file <- read.csv(ipfile1$datapath, header = FALSE, stringsAsFactors=FALSE)
            new.break.file <- as.numeric(break.file[1,])
            ipfile2 <- input$file2
            labels.file <- read.csv(ipfile2$datapath, header = FALSE, stringsAsFactors=FALSE)
            new.labels.file <- as.numeric(labels.file[1,])
            plot_HorF(input$allele, y = input$choose, breaks = new.break.file, label =  new.labels.file)}
        }
        
        else{
          if(is.null(input$file1)){new.break.file = NULL
          if(is.null(input$file2)){new.labels.file = NULL}
          plot_HorF(input$allele, y = input$choose)}  # file1 should be a .txt file
          else{
            ipfile1 <- input$file1
            break.file <- read.csv(ipfile1$datapath, header = FALSE, stringsAsFactors=FALSE)
            new.break.file <- as.numeric(break.file[1,])
            ipfile2 <- input$file2
            labels.file <- read.csv(ipfile2$datapath, header = FALSE, stringsAsFactors=FALSE)
            new.labels.file <- as.numeric(labels.file[1,])
            plot_HorF(input$allele, y = input$choose, breaks = new.break.file, label =  new.labels.file)}
        }
        
      })
      
      output$dynamic_value1 <- renderPrint({
        if(input$choose == "Hadm"){
          if(is.null(input$file1)){new.break.file = NULL }  # file1 should be a .txt file
          else{
            ipfile1 <- input$file1
            break.file <- read.csv(ipfile1$datapath, header = FALSE, stringsAsFactors=FALSE)
            new.break.file <- as.numeric(break.file[1,])}
        }
        if(input$choose == "Fst"){
          if(is.null(input$file1)){new.break.file = NULL}  # file1 should be a .txt file
          else{
            ipfile1 <- input$file1
            break.file <- read.csv(ipfile1$datapath, header = FALSE, stringsAsFactors=FALSE)
            new.break.file <- as.numeric(break.file[1,])}
        }
        str(new.break.file)
      })
      
      output$dynamic_value2 <- renderPrint({
        if(input$choose == "Hadm"){
          if(is.null(input$file2)){new.labels.file = NULL}  # file1 should be a .txt file
          else{
            ipfile2 <- input$file2
            labels.file <- read.csv(ipfile2$datapath, header = FALSE, stringsAsFactors=FALSE)
            new.labels.file <- as.numeric(labels.file[1,])}
        }
        if(input$choose == "Fst"){
          if(is.null(input$file2)){new.labels.file = NULL}  # file1 should be a .txt file
          else{
            ipfile2 <- input$file2
            labels.file <- read.csv(ipfile2$datapath, header = FALSE, stringsAsFactors=FALSE)
            new.labels.file <- as.numeric(labels.file[1,])}
        }
        str(new.labels.file)
      })
      
    })
}

shinyApp(ui, server)