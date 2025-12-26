# Load in libraries
library(shiny)
library(tidyverse)
library(rlang)

# Read in defined functions
source("cov_matrix.R")
source("power_calc.R")
source("plot_power.R")
source("workflow.R")

ui <- fluidPage(
  
  titlePanel("MMRM Power / Sample Size Calculator"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("cov_type", "Covariance type", 
                  choices = c("AR", "CS"), 
                  selected = "AR"),
      numericInput("rho", "Correlation parameter", 
                   value = 0.8, 
                   min = 0, max = 1, step = 0.05),
      textInput("times", "Visit times", 
                value = "8, 18, 28, 40, 52"),
      textInput("stds", "Standard deviations", 
                value = "175, 200, 225, 250, 275"),
      textInput("dropout_list", "Dropout rate(s)", 
                value = "0, 0.25, 0.3"),
      textInput("margin_list", "Margin(s)", 
                value = "-40, -45, -50"),
      textInput("diff_list", "Treatment effect(s)", 
                value = "20"),
      textInput("power_list", "Target power(s) [leave blank to calculate]", 
                value = "0.8, 0.9"),
      textInput("ntot_list", "Total sample size(s) [leave blank to calculate]", 
                value = ""),
      selectInput("type", "Trial type", 
                  choices = c("SUP", "NI", "EQUI"), 
                  selected = "NI"),
      numericInput("q", "Number of model covariates", 
                   value = 3, min = 0),
      numericInput("alpha", "Alpha", 
                   value = 0.05),
      numericInput("gamma0", "Proportion in control group", 
                   value = 0.5),
      selectInput("plot_yn", "Generate plots?", 
                  choices = c("Y", "N"), 
                  selected = "N"),
      textInput("xlabel", "X-axis label (optional)", 
                value = ""),
      textInput("ylabel", "Y-axis label (optional)", 
                value = ""),
      textInput("title", "Plot title (optional)", 
                value = ""),
      actionButton("run", "Run")
    ),
    
    mainPanel(
      verbatimTextOutput("matrix_out"),
      verbatimTextOutput("res_table"),
      uiOutput("plot_output")
    )
  )
)

server <- function(input, output, session) {
  
  parse_numeric_list <- function(x) {
    x <- gsub(" ", "", x)
    if (x == "") return(NA)
    if (tolower(x) == "na") return(NA)
    as.numeric(unlist(strsplit(x, ",")))
  }
  
  run_results <- eventReactive(input$run, {
    times <- parse_numeric_list(input$times)
    stds <- parse_numeric_list(input$stds)
    dropout_list <- parse_numeric_list(input$dropout_list)
    margin_list <- parse_numeric_list(input$margin_list)
    diff_list <- parse_numeric_list(input$diff_list)
    power_list <- parse_numeric_list(input$power_list)
    ntot_list <- parse_numeric_list(input$ntot_list)
    
    word <- NULL
    if (all(is.na(ntot_list))) {
      word <- "sample size"
    } else if (all(is.na(power_list))) {
      word <- "power"
    } else {
      word <- ""
    }
    
    out <- run_workflow(
      times = times,
      stds = stds,
      cov_type = input$cov_type,
      rho = input$rho,
      dropout_list = dropout_list,
      margin_list = margin_list,
      diff_list = diff_list,
      power_list = power_list,
      ntot_list = ntot_list,
      type = tolower(input$type),
      q = input$q,
      alpha = input$alpha,
      gamma0 = input$gamma0,
      plot_yn = input$plot_yn,
      title = ifelse(nzchar(input$title), input$title, NA),
      xlabel = ifelse(nzchar(input$xlabel), input$xlabel, NA),
      ylabel = ifelse(nzchar(input$ylabel), input$ylabel, NA)
    )
    
    out$word <- word
    
    return(out)
    
  })
  
  output$matrix_out <- renderPrint({
    req(run_results())
    cat(paste0(toupper(input$cov_type), "(", input$rho, ") covariance matrix:\n"))
    print(run_results()$covmat)
  })
  
  output$res_table <- renderPrint({
    
    req(run_results())
    tbl <- run_results()$table[,-6]
    
    if ("ntot" %in% names(tbl)) {
      tbl$ntot <- round(tbl$ntot)
    } else if (is.numeric(tbl[[1]])) {
      tbl[[1]] <- round(tbl[[1]])
    }
    
    word <- run_results()$word
    cat(paste0(toupper(input$type), " ", word, " results:\n"))
    print(tbl, row.names = FALSE)
    
  })
  
  output$plot_output <- renderUI({
    req(run_results())
    if (toupper(input$plot_yn) == "Y" && !is.null(run_results()$plots)) {
      valid_plots <- run_results()$plots[!sapply(run_results()$plots, is.null)]
      plot_output_list <- lapply(seq_along(valid_plots), function(i) {
        plotname <- paste0("plot_", i)
        plotOutput(plotname, height = "400px", width = "600px")
      })
      do.call(tagList, plot_output_list)
    }
  })

  observe({
    req(run_results())
    if (toupper(input$plot_yn) == "Y" && !is.null(run_results()$plots)) {
      valid_plots <- run_results()$plots[!sapply(run_results()$plots, is.null)]
      
      lapply(seq_along(valid_plots), function(i) {
        local({
          my_i <- i
          plotname <- paste0("plot_", my_i)
          output[[plotname]] <- renderPlot({
            valid_plots[[my_i]]
          })
        })
      })
    }
  })
  
}

# Launch app
shinyApp(ui = ui, server = server)