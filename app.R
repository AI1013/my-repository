library(shiny)
source("appfunction.R", local = TRUE)
ui <- fluidPage(
  sidebarPanel(textInput("allele", "Enterloci"),
  checkboxGroupInput("choose", "Please choose", choices = c("Hadm", "Fst"), selected = NULL),
  fileInput("file1", "Set breaks"),
  fileInput("file2", "Set labels"),
  actionButton("BUTTON1","Press the button")),
  
  
  mainPanel(
  plotOutput("Plot", width = "100%", height = "400px"))
)

server <- function(input, output, session) {
  
  observeEvent(
  input$BUTTON1, 
  {
  output$Plot <- renderPlot({
      
      if(is.null(input$file1)){return(NULL)} # file1 should be a .txt file
      ipfile1 <- input$file1
      break.file <- read.csv(ipfile1$datapath, header = FALSE, stringsAsFactors=FALSE)
      new.break.file <- as.numeric(break.file[1,])
    
  
      if(is.null(input$file2)){return(NULL)} # file2 should be a .txt file
      ipfile2 <- input$file2
      labels.file <- read.csv(ipfile2$datapath, header = FALSE, stringsAsFactors=FALSE)
      new.labels.file <- as.numeric(labels.file[1,]) 
    
    
    plot_HorF(input$allele, y = input$choose, breaks = new.break.file, label =  new.labels.file)

  })})
}

shinyApp(ui, server)