library(shiny)
library(markdown)
library(knitr)

fluidPage(
  headerPanel("Detection of Significance Reversers in Linear Regression"),
  fluidRow(column(2, 
    ## include text in readme.md
    includeMarkdown("readme.md"),
    ## example button
    includeMarkdown("Example.md"),
    p("(Figure of a 2015 PNAS paper)"),
    actionButton("run.example", "Example"),
    br(), br(),
    ## file input
    includeMarkdown("Import.md"),
    ## header selection
    checkboxInput("header", "Data has a Header row", TRUE),
    ## csv type selection
    radioButtons("dec.type", "Type of decimal dot",
                 choices = c("dot (.)" = "dec1", "comma (,)" = "dec2"),
                 selected = "dec1", inline = TRUE),
    radioButtons("sep.type", "Type of separator",
                 choices = c("comma (,)" = "sep1", "space (\\s)" = "sep2",
                             "tab (\\t)" = "sep3", "semicolon (;)" = "sep4"), 
                 selected = "sep1", inline = TRUE),
    ## file input
    fileInput("input.file", "Choose .csv file",
              accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
    ## significance border
    textInput("pval", "P-value significance border:", 0.05),
    ## download plot, stats and Cook's D
    includeMarkdown("Results.md"),
    downloadButton("plot.download", "Download Plots"),   
    downloadButton("stat.download", "Download Stats")
  ),
  ## main panel with "dynamic.tabset" output
  column(10, 
      uiOutput("dynamic.tabset")
  )
)
)