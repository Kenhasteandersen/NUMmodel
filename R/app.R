#
# Install app:
#  ssh ken@oceanlife.dtuaqua.dk
#  cd FirstPrinciplesPlankton
#  update the git (git pull)
#  sudo cp -R ~/FirstPrinciplesPlankton/*  /srv/shiny-server/Plankton2
#  If using new packages install them by running R as root (sudo su; R; install.packages("XXX))
#  sudo systemctl restart shiny-server
# 

library(shiny)
options(shiny.sanitize.errors = FALSE)

source("modelChemostat.R")

bUseC = FALSE
bUseF = TRUE

# ===================================================
# Define UI for application
# ===================================================

uiChemostat <- fluidPage(
  tags$head(
     # Add google analytics tracking:
     includeHTML(("analytics.html")),
     # Make rules widers:
     tags$style(HTML("hr {border-top: 1px solid #000000;}"))
  )
  ,
  h1('Size-based plankton simulator'),
  p('Simulate a plankton ecosystem in the upper part of a watercolumn. 
   Cell size is the only trait characterizing each plankton group.
    All groups are able to perform photoharvesting, taking up dissolve nutrients and carbon, and do phagotrophy.
    The trophic strategy (the background colours) is an emergent property.'),
  p('Documentation in:',
    a("Andersen and Visser (2023)",href="https://www.biorxiv.org/content/10.1101/2022.05.16.492092v2"),'.'),
  p('Revised version 1.1, January 2023.')
  ,
  # Sidebar with a slider inputs
  sidebarLayout(
    sidebarPanel(
      selectInput("setup", "Setup:",
                  c("Unicellular (simple)" = "GeneralistsSimpleOnly",
                    "Unicellular" = "GeneralistsOnly"))
                    #"Uni+multicellular" = "Generic"))
      ,
      sliderInput("L",
                  "Light (PAR; uE/m2/s)",
                  min = 0,
                  max = 300,
                  step=1,
                  value = 60)
      ,
      sliderInput("d10",
                  "log10(Mixing rate) (1/day)",
                  min = -4,
                  max = 0,
                  step = 0.1,
                  value = log10(parametersChemostat()$d))
      ,
      sliderInput("T",
                  "Temperature",
                  min = 0,
                  max = 25,
                  step = 0.5,
                  value = 10)
      ,
      #sliderInput("latitude",
      #            "Latitude (degrees; experimental)",
      #            min = 0,
      #            max = 90,
      #            step=1,
      #            value = 0)
      #,
      hr()
      ,
      #sliderInput("N0",
      #            "Deep nutrients (mugN/l)",
      #            min = 0,
      #            max = 200,
      #            step = 0.1,
      #            value = 14)
      sliderInput("mortHTL",
                  "Mortality by higher trophic levels on the largest sizes (1/day)",
                  min = 0,
                  max = 0.5,
                  step = 0.01,
                  value = parametersChemostat()$mortHTL)
      ,
      # sliderInput("mort2",
      #             "Quadratic mortality coef. (virulysis; 1/day/mugC)",
      #             min = 0,
      #             max = 0.01,
      #             step = 0.0001,
      #             value = 0.004)
      # ,
      # sliderInput("epsilon_r",
      #             "Remineralisation of virulysis and HTL losses",
      #             min = 0,
      #             max = 1,
      #             value = 0,
      #             step = 0.1)
      # ,
      sliderInput("M",
                  "Thickness of productive layer (m)",
                  min = 0,
                  max = 100,
                  step = 1,
                  value = parametersChemostat()$M)
      
    ),
    #
    # Show plots:
    #
    mainPanel(
      #
      # Main output panel
      #
      tabsetPanel(
        tabPanel(
          "Main output",
          #conditionalPanel(
          #  condition = "input.latitude>0",
          #  sliderInput("t",
          #              "Time (day)",
          #              min=0, 
          #              max=365,
          #              step=1,
          #              value=120,
          #              width="600px"),
          #  plotOutput("plotSeasonal", width="600px", height="100px")  
          #),
          plotOutput("plotSpectrum", width="600px", height="300px"),
          plotOutput("plotRates", width="600px", height="300px"),
          plotOutput("plotLeaks", width="600px", height="200px"),
          plotOutput("plotDeltas", width="600px", height="200px")
        )
        ,
        tabPanel(
          "Ecosystem functions",
          plotOutput("plotFunction", width="600px")
        )
        ,
        tabPanel(
          "Dynamics",
          conditionalPanel(
            condition="input.latitude>0",
            plotOutput("plotSeasonalTimeline", width="600px")
          ),
          plotOutput("plotTime", width="600px")#,
          #plotOutput("plotComplexRates", width="700px")
        )
      )
    )))
#
# Define server logic
#
serverChemostat <- function(input, output) {
  #
  # Simulate the system when a parameter is changed:
  #
  sim <- eventReactive({
    input$setup
    input$L
    input$latitude
    input$d10
    input$T
    input$M
    #input$mort2
    input$mortHTL
    #input$epsilon_r
  },
  {
    # get all base parameters
    p <- parametersChemostat( parameters(n=20) ) 
    # Update parameters with user input:
    p$L = input$L
    p$latitude = 0#input$latitude
    p$d = 10^input$d10
    p$T = input$T
    #p$mort2 = input$mort2/log(p$Delta)
    p$mortHTL = input$mortHTL
    p$M = input$M
    #p$remin = input$epsilon_r
    
    #p$tEnd = 1000
    
    #if (p$latitude>0)
    #  p$tEnd = 2*365

    # Simulate
    return(simulateChemostat(p, bUseC, bUseF, input$setup))   
  })
  #
  # Plots:
  #
  output$plotSpectrum <- renderPlot(plotSpectrum(sim(), input$t))
  output$plotRates <- renderPlot(plotRates(sim(), t=input$t))
  output$plotLeaks = renderPlot(plotLeaks(sim(), input$t))
  output$plotDeltas = renderPlot(plotDeltas(sim()))
  output$plotComplexRates <- renderPlot(plotComplexRates(sim(), input$t))
  output$plotTime <- renderPlot(plotTimeline(sim(), input$t))
  output$plotSeasonalTimeline <- renderPlot(plotSeasonalTimeline(sim()))
  output$plotFunction <- renderPlot(plotFunctionsChemostat(sim()))
  output$plotSeasonal <- renderPlot(plotSeasonal(p, input$t))
}

# Compile C code:
#if (system2("g++", "-O3 -shared ../Cpp/model.cpp -o ../Cpp/model.so")==0)  {
#  bUseC = TRUE
#  cat("Using C engine\n")
#}

# Run the application 
shinyApp(ui = uiChemostat, server = serverChemostat)

