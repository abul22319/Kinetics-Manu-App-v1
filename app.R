options(shiny.maxRequestSize = 30*1024^2, scipen = 999,
        shiny.launch.browser = .rs.invokeShinyWindowExternal)

source("Monoexponential.R")

packages = c(
  "shiny",
  "shinydashboard",
  "shinythemes",
  "shinyBS",
  "tidyverse",
  "ggplot2",
  "lme4",
  "readxl",
  "signal",
  "zoo",
  "writexl",
  "dplyr",
  "bslib",
  "plotly",
  "shinyjs",
  "minpack.lm"
)

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

library(shiny)


# =============================================================================
# UI
# =============================================================================

ui <- page_sidebar(
  theme = bs_theme(preset = "simplex"),
  title = "Exponential Model Fitting",
  
  sidebar = sidebar(
    width = 300,
    
    fileInput("file",
              label = "Select File"),
    
    radioButtons("direction",
                 label   = "Choose Model Direction:",
                 choices = list("Rise"  = 1,
                                "Decay" = 2),
                 selected = 1),
    
    radioButtons("model_type",
                 label   = "Choose Model Type:",
                 choices = list("1-Component (Mono-exponential)" = "mono",
                                "2-Component (Bi-exponential)"   = "bi"),
                 selected = "mono"),
    
    hr(),
    
    sliderInput("fc",
                label   = "Low-pass filter cutoff (Hz):",
                min     = 0.01,
                max     = 0.16,
                value   = 0.15,
                step    = 0.01),
    
    helpText("Sampling rate is 0.333 Hz → Nyquist = 0.167 Hz.
              Lower values apply more smoothing before model fitting."),
    
    hr(),
    
    # Dynamic info box explaining the selected model
    uiOutput("model_info")
  ),
  
  navset_bar(
    
    # ---- Absolute FBF ----
    nav_panel("Absolute FBF",
              fluidRow(
                column(width = 6,
                       card(full_screen = TRUE,
                            card_header("Model Fit"),
                            plotOutput("Model.Plot")),
                       card(card_header("Parameters"),
                            tableOutput("Parameters"))
                ),
                column(width = 6,
                       card(full_screen = TRUE,
                            card_header("Residuals"),
                            plotOutput("Line.Plot")),
                       card(card_header("Correlation Results"),
                            tableOutput("Cor.Result"))
                )
              )
    ),
    
    # ---- Relative FBF ----
    nav_panel("Relative FBF",
              fluidRow(
                column(width = 6,
                       card(full_screen = TRUE,
                            card_header("Model Fit"),
                            plotOutput("Model.Plot_FBF_rel")),
                       card(card_header("Parameters"),
                            tableOutput("Parameters_FBF_rel"))
                ),
                column(width = 6,
                       card(full_screen = TRUE,
                            card_header("Residuals"),
                            plotOutput("Line.Plot_FBF_rel")),
                       card(card_header("Correlation Results"),
                            tableOutput("Cor.Result_FBF_rel"))
                )
              )
    ),
    
    # ---- Absolute FVC ----
    nav_panel("Absolute FVC",
              fluidRow(
                column(width = 6,
                       card(full_screen = TRUE,
                            card_header("Model Fit"),
                            plotOutput("Model.Plot_FVC")),
                       card(card_header("Parameters"),
                            tableOutput("Parameters_FVC"))
                ),
                column(width = 6,
                       card(full_screen = TRUE,
                            card_header("Residuals"),
                            plotOutput("Line.Plot_FVC")),
                       card(card_header("Correlation Results"),
                            tableOutput("Cor.Result_FVC"))
                )
              )
    ),
    
    # ---- Relative FVC ----
    nav_panel("Relative FVC",
              fluidRow(
                column(width = 6,
                       card(full_screen = TRUE,
                            card_header("Model Fit"),
                            plotOutput("Model.Plot_FVC_rel")),
                       card(card_header("Parameters"),
                            tableOutput("Parameters_FVC_rel"))
                ),
                column(width = 6,
                       card(full_screen = TRUE,
                            card_header("Residuals"),
                            plotOutput("Line.Plot_FVC_rel")),
                       card(card_header("Correlation Results"),
                            tableOutput("Cor.Result_FVC_rel"))
                )
              )
    )
  )
)


# =============================================================================
# SERVER
# =============================================================================

server <- function(input, output, session) {
  
  session$onSessionEnded(function() {
    stopApp()
  })
  
  # ---- Dynamic model info box in sidebar ----
  output$model_info <- renderUI({
    if (input$model_type == "mono") {
      tags$div(
        style = "background:#f0f4ff; border-left:4px solid #4a6fa5;
                 padding:10px; border-radius:4px; font-size:0.85em;",
        tags$b("Mono-exponential formula:"),
        tags$br(),
        tags$code("Y ~ B * (1 - exp(-(t - TD) / tau))"),
        tags$br(), tags$br(),
        tags$b("Parameters:"),
        tags$ul(
          tags$li("B — Amplitude"),
          tags$li("tau — Time constant"),
          tags$li("TD — Time delay"),
          tags$li("MRT — Mean response time (tau + TD)")
        )
      )
    } else {
      tags$div(
        style = "background:#fff4f0; border-left:4px solid #a54a4a;
                 padding:10px; border-radius:4px; font-size:0.85em;",
        tags$b("Bi-exponential formula:"),
        tags$br(),
        tags$code("Y ~ B1*(1-exp(-(t-TD1)/tau1)) + B2*(1-exp(-(t-TD2)/tau2))"),
        tags$br(), tags$br(),
        tags$b("Parameters:"),
        tags$ul(
          tags$li("B1 / B2 — Component amplitudes"),
          tags$li("tau1 / tau2 — Time constants"),
          tags$li("TD1 / TD2 — Time delays"),
          tags$li("MRT1 / MRT2 — Mean response times")
        )
      )
    }
  })
  
  
  # ---- Core reactive: run selected model ----
  Work2BDone <- reactive({
    
    req(input$file$datapath, input$direction, input$model_type, input$fc)
    
    withProgress(message = "Fitting model, please wait...", value = 0.5, {
      
      if (input$model_type == "mono") {
        MonoExp(input$file$datapath, input$direction, fc = input$fc)
      } else {
        BiExp(input$file$datapath, input$direction, fc = input$fc)
      }
    })
  })
  
  
  # ---- Absolute FBF outputs ----
  output$Model.Plot <- renderPlot({
    req(Work2BDone())
    Work2BDone()$Exp.Model_FBF
  })
  
  output$Line.Plot <- renderPlot({
    req(Work2BDone())
    Work2BDone()$RefLine.Model_FBF
  })
  
  output$Parameters <- renderTable({
    req(Work2BDone())
    Work2BDone()$Parameters_FBF
  }, digits = 1)
  
  output$Cor.Result <- renderTable({
    req(Work2BDone())
    Work2BDone()$Cor.Result_FBF
  }, digits = 1)
  
  
  # ---- Relative FBF outputs ----
  output$Model.Plot_FBF_rel <- renderPlot({
    req(Work2BDone())
    Work2BDone()$Exp.Model_FBF_rel
  })
  
  output$Line.Plot_FBF_rel <- renderPlot({
    req(Work2BDone())
    Work2BDone()$RefLine.Model_FBF_rel
  })
  
  output$Parameters_FBF_rel <- renderTable({
    req(Work2BDone())
    Work2BDone()$Parameters_FBF_rel
  }, digits = 1)
  
  output$Cor.Result_FBF_rel <- renderTable({
    req(Work2BDone())
    Work2BDone()$Cor.Result_FBF_rel
  }, digits = 1)
  
  
  # ---- Absolute FVC outputs ----
  output$Model.Plot_FVC <- renderPlot({
    req(Work2BDone())
    Work2BDone()$Exp.Model_FVC
  })
  
  output$Line.Plot_FVC <- renderPlot({
    req(Work2BDone())
    Work2BDone()$RefLine.Model_FVC
  })
  
  output$Parameters_FVC <- renderTable({
    req(Work2BDone())
    Work2BDone()$Parameters_FVC
  }, digits = 1)
  
  output$Cor.Result_FVC <- renderTable({
    req(Work2BDone())
    Work2BDone()$Cor.Result_FVC
  }, digits = 1)
  
  
  # ---- Relative FVC outputs ----
  output$Model.Plot_FVC_rel <- renderPlot({
    req(Work2BDone())
    Work2BDone()$Exp.Model_FVC_rel
  })
  
  output$Line.Plot_FVC_rel <- renderPlot({
    req(Work2BDone())
    Work2BDone()$RefLine.Model_FVC_rel
  })
  
  output$Parameters_FVC_rel <- renderTable({
    req(Work2BDone())
    Work2BDone()$Parameters_FVC_rel
  }, digits = 1)
  
  output$Cor.Result_FVC_rel <- renderTable({
    req(Work2BDone())
    Work2BDone()$Cor.Result_FVC_rel
  }, digits = 1)
  
}


# =============================================================================
# RUN
# =============================================================================

shinyApp(ui = ui, server = server)