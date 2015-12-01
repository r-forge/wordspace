##
##  Shiny semantic similarity GUI -- GUI definition
##

shinyUI(pageWithSidebar(
  headerPanel("Semantic similarity explorer"),

  sidebarPanel(
    textInput("target", "Find nearest neighbours of:", value="", width="100%"),
    selectInput("candidates", label=NULL, choices=NULL, size=8, selectize=FALSE, width="100%"),
    sliderInput("NN", "Number of neighbours", min=5, max=50, value=20, step=1, width="100%"),
    tags$div(style = "margin-top: 1em;",
      tags$span(style = "font-weight: bold;", "Instructions:"),
      tags$ul(
        tags$li("Enter complete target word at the top, or enter prefix and select from list of completions below."),
        tags$li("Append", tags$code("%"), "for explicit prefix search."),
        tags$li("Click on node in neighbour graph to select word as new target.")
        )
      )
    ),

  mainPanel(
    tabsetPanel(
      tabPanel("Graph",
        fluidRow(
          column(3, checkboxInput("isoMDS", "non-metric MDS", value=TRUE)),
          column(3, checkboxInput("edges", "show edges", value=TRUE)),
          column(6, sliderInput("edgeWidth", NULL, min=1, max=20, value=6, step=.1, ticks=FALSE))
          ),
        plotOutput("NNplot", width="100%", height="600px", click="mapClick")),
      tabPanel("List",
               htmlOutput("NNlist")),
      tabPanel("Debug", 
               textOutput("debug", container = function (...) pre(..., style="font-size: 75%;")))
      )
    )
  ))

      

