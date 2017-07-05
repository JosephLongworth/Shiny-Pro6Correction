################################################################################################
################################################################################################
################################################################################################
# This is the ui script for the Proline 6 correction Shiny app developed by Joseph Longworth

################################################################################################
################################################################################################
################################################################################################
library(shiny)
library(stringr)
shinyUI(fluidPage(

  titlePanel("Proline 6 correction by JLongworth (j.longworth@sheffield.ac.uk)"),
  sidebarLayout(
    sidebarPanel(
#      textInput("Name", "Name:", "Submitter name"),
#      textInput("run_name", "Name for search:", "HELA 200g"),
      fileInput('file1', 'select protein groups file',
                 accept=c('text/txt', 
                          'text/comma-separated-values,text/plain', 
                          '.txt')),
      fileInput('file2', 'select evidence file',
                accept=c('text/txt', 
                         'text/comma-separated-values,text/plain', 
                         '.txt')),

      tags$hr()
            ),
    mainPanel(
                        plotOutput("plot"),
                        uiOutput("download_button")
    
    )
                        )
))