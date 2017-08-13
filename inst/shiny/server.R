library(shiny)
library(markdown)
library(knitr)
#source("InfluenceApp.R")

## server for the Shiny app
shinyServer(function(input, output, session) {
  
  ## input for x/y data or 'Example'
  imported.data <- reactive({
    if(input[["run.example"]] >= 1) impdat <- read.csv("XY-data.csv")
   
    if (!is.null(input[["input.file"]])) {  
      header <- input[["header"]]
      dec <- switch(input[["dec.type"]], dec1 = ".", dec2 = ",")
      sep <- switch(input[["sep.type"]], sep1 = ",", sep2 = " ", sep3 = "\t", sep4 = ";")
      
      impdat <- read.csv(input[["input.file"]][["datapath"]], header = header, 
                      sep = sep, dec = dec)
      
      if(input[["header"]] == FALSE) {
        if (ncol(impdat) == 2) colnames(impdat) <- c("x", "y")
        if (ncol(impdat) == 3) colnames(impdat) <- c("x", "y", "w")
      }
    }
    if (NCOL(impdat) == 1) return(1) else return(impdat)
  })
  
  ## dynamic tabset before and after data input
  output[["dynamic.tabset"]] <- renderUI({
    if(input[["run.example"]] == 0 && is.null(input[["input.file"]])) {
      tabPanel("No input detected", h3("No input detected! Run 'Example' or import data."))
    } else {
      tabsetPanel(id = "tabs", type = "pills",
        tabPanel("Statistics", 
                 p("\b"),
                 h4("Message:"),
                 h5("'Error: incorrect number of dimensions' indicates wrong import parameters."),
                 verbatimTextOutput("message.output"),
                 p("\b"),
                 h4("Summary of full model (all data points):"),
                 renderPrint({summary(res.infl()$origModel)}),
                 p("\b"),
                 h4("Statistical Analysis for each Leave-One-Out model. P-value reversers are marked. For abbreviations, see below Table."),
                 output$DT <- DT::renderDataTable(signif(res.infl()[["infl"]], 5), 
                   options = list(paging = FALSE, searching = FALSE),
                   rownames = FALSE),
                 h4("Legend:"),
                 h5(strong("dSlope:"), "      difference to original slope."),
                 h5(strong("dInter:"), "      difference to original intercept."),
                 h5(strong("dfb.Inter:"), "   dfbeta of intercept."),
                 h5(strong("dfb.Slope:"), "   dfbeta of slope."),
                 h5(strong("dffit:"), "       dffits."),
                 h5(strong("cov.r:"), "       covariance ratio."),
                 h5(strong("cook.d:"), "      Cook's distance."),
                 h5(strong("dfb.Slope:"), "   dfbeta of slope."),
                 h5(strong("hat:"), "         diagonal of hat matrix."),
                 h5(strong("looP:"), "        leave-one-out p-value."),
                 h5(strong("dP:"), "          difference to full model p-value."),
                 h5(strong("SE:"), "          standard error of slope."),
                 h5(strong("dSE:"), "         difference to full model standard error of slope."),
                 h5(strong("Rsq:"), "         R-square.")),
        tabPanel("Regression Plot", 
                 p("\b"),
                 h4("Linear Regression Plot with Influential Data."),
                 h4("Significance reversers are depicted with larger red points!"),
                 h4("Mouse over data point to see its x/y value."),
                 p("\b"),
                 verbatimTextOutput("hover_info"),
                 plotOutput("lmPlot.output", hover = hoverOpts(id = "plot_hover"))),
        tabPanel("P-Value Plot", 
                 p("\b"),
                 h4("P-value Plot for each Leave-One-Out Data Point."),
                 h4("Significance reversers will cross the red line!"),
                 p("\b"),
                 plotOutput("pvalPlot.output")),
        tabPanel("Slope vs. SE Plot", 
                 p("\b"),
                 h4("Slope versus S.E. Plot for each Leave-One-Out Data Point."),
                 h4("Leave-One-Out of points on the right side of the t-border result in a significant model, 
                    while points on the left side result in an insignificant model!"),
                 h4("Dashed lines represent the slope & s.e.(slope) of the full data model."),
                 plotOutput("slsePlot.output")),
        tabPanel("Influence Plot", 
                 p("\b"),
                 h4("Influence Measures versus delta-P values Plot for each Leave-One-Out Data Point."),
                 h4("Significance reversers are depicted as red dots!"),
                 h4("Often, one can observe significance reversal for data points with Cook's distance/leverage/dfbeta(slope)/dffits/covratio below the commonly applied cut-off values."),
                 p("\b"),
                 plotOutput("inflPlot.output"))
      )
    }
  })
  
  ## Run lmInfl function
  res.infl <- reactive({   
    DAT <- imported.data() 
    if (NCOL(DAT) == 3) wts <- DAT[, 3] else wts <- NULL
    MODEL <- lm(DAT[, 2] ~ DAT[, 1], weights = wts)
    colnames(MODEL$model) <- colnames(DAT)[2:1]
    res <- lmInfl(model = MODEL, alpha = as.numeric(input[["pval"]]), verbose = FALSE) 
    return(res)
  })
  
  ## lmPlot output
  output[["lmPlot.output"]] <- renderPlot({
    lmPlot(res.infl())
  }, width = 800, height = 600)
  
  ## hover over lmPlot
  output[["hover_info"]] <- renderPrint({
    cat("Data under cursor:\n")
    NP <- nearPoints(res.infl()$infl[, c(1:5, 14)], input$plot_hover, 
                    xvar = colnames(res.infl()$infl)[2],
                    yvar = colnames(res.infl()$infl)[3], maxpoints = 1)
    if (nrow(NP) == 0) cat("No data point.\n   ") else print(NP)
  })
  
  ## pvalPlot output
  output[["pvalPlot.output"]] <- renderPlot({
    pvalPlot(res.infl())
  }, width = 800, height = 600)
  
  ## slsePlot output
  output[["slsePlot.output"]] <- renderPlot({
    slsePlot(res.infl())
  }, width = 800, height = 600)
  
  ## inflPlot output
  output[["inflPlot.output"]] <- renderPlot({
    inflPlot(res.infl())
  }, width = 1000, height = 600)
  
  ## Message output
  output[["message.output"]] <- renderPrint({
    RF <- res.infl()
    if ((RF$origP <= RF$alpha) & length(RF$sel) > 0) {
      cat("Significant model with p =", summary(RF$origModel)$coefficients[2, 4],
          "has been changed to insignificant by removing one of the following data points (in order of strength):\n\n")
      PVALS <- RF$infl[RF$sel, "looP"]
      for (i in 1:length(PVALS)) cat("#", RF$sel[i], " => ", PVALS[i], "\n", sep = "") 
    }
    if ((RF$origP > RF$alpha) & length(RF$sel) > 0) {
      cat("Insignificant model with p =", summary(RF$origModel)$coefficients[2, 4],
          "has been changed to significant by removing one of the following data points (in order of strength):\n\n")
      PVALS <- RF$infl[RF$sel, "looP"]
      for (i in 1:length(PVALS)) cat("#", RF$sel[i], " => ", PVALS[i], "\n", sep = "") 
    }
    if ((RF$origP <= RF$alpha) & length(RF$sel) == 0) {
      cat("Significant model with p =", summary(RF$origModel)$coefficients[2, 4],
          "has no significance reversers.")
    }
    if ((RF$origP > RF$alpha) & length(RF$sel) == 0) {
      cat("Insignificant model with p =", summary(RF$origModel)$coefficients[2, 4],
          "has no significance reversers.")
    }
  })
  
  ## Plots download
  output[["plot.download"]] <- downloadHandler(
    filename = "Plots.pdf",
    content = function(file) {
      pdf(file)
      lmPlot(res.infl())
      pvalPlot(res.infl())
      slsePlot(res.infl())
      inflPlot(res.infl())
      dev.off()    
  })  
  
  output[["stat.download"]] <- downloadHandler(
    filename = "Stats.txt",
    content = function(file) {
      sink(file)
      #### Full model with all data points:
      print("Full model with all data points:") 
      print(summary(res.infl()$origModel))
      #### Stats
      print("Influence Measures for all LOO-models:")
      print(res.infl()$infl)
      sink()
    })  
  
})
