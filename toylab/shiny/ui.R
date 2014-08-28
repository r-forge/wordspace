##
##  Shiny semantic similarity GUI -- GUI definition
##

shinyUI(pageWithSidebar(
  headerPanel("Semantic similarity explorer"),

  sidebarPanel(
    ## selectizeInput("target", "Find nearest neighbours of:", choices=sort(rownames(M)), width="100%",
    ##                options=list(placeholder="please type a word", onInitialize = I('function() { this.setValue(""); }'))),
    selectizeInput("target", "Find nearest neighbours of:", choices=NULL, multiple=FALSE, width="100%",
                   options=list(placeholder="enter search term", maxOptions=5)),
    sliderInput("NN", "number of neighbours", min=5, max=50, value=20, step=1, width="100%")
    ),

  mainPanel(
    tabsetPanel(
      tabPanel("Graph",
               tags$table(tags$tr(
                 tags$td(checkboxInput("isoMDS", "non-metric MDS", value=TRUE)),
                 tags$td(checkboxInput("edges", "show edges", value=TRUE), style="padding-left: 2em;"),
                 tags$td(sliderInput("edgeWidth", NULL, min=1, max=20, value=6, step=.1, ticks=FALSE), style="padding-left: 2em;")
                 )),
               plotOutput("NNplot", width="600px", height="600px", clickId="mapClick")),
      tabPanel("List",
               htmlOutput("NNlist")),
      tabPanel("Debug", 
               textOutput("debug", container = function (...) pre(..., style="font-size: 75%;")))
      )
    )
  ))

      

