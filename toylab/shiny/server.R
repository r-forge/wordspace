##
##  Shiny semantic similarity GUI -- Server implementation
##

shinyServer(function(input, output, session) {
  terms.sorted <- sort(rownames(M))
  TargetInfo <- data.frame(label=terms.sorted, value=terms.sorted)
  updateSelectizeInput(session, "target", choices = TargetInfo, server = TRUE)

  ENV <- reactiveValues(coord=NULL)
  
  NN <- reactive({
    target <- input$target
    ok <- target %in% rownames(M)
    if (ok) {
      nearest.neighbours(M, target, n=input$NN, dist.matrix=TRUE)
    } else {
      NULL
    } 
  })
  
  output$NNplot <- renderPlot({
    nn.dist <- NN()
    par(mar=c(0,0,0,0))
    if (is.null(nn.dist)) {
      plot(0, 0, type="n", xlim=c(0, 2), ylim=c(0, 2), xaxt="n", yaxt="n", xlab="", ylab="")
      text(1, 1, sprintf("%s not found in DSM", input$target))
    } else {
      ENV$coord <- plot(nn.dist, expand=.1, cex=1.2, method=if (input$isoMDS) "isomds" else "sammon", show.edges=input$edges, edges.lwd=input$edgeWidth)
    }
  })

  output$NNlist <- renderUI({
    nn.dist <- NN()
    if (is.null(nn.dist)) {
      span(tags$em(input$target), "not found in DSM", style="color: red;")
    } else {
      nn.terms <- rownames(nn.dist)[-1]
      list.items <- lapply(nn.terms, function (.t) {
        tags$li(.t, tags$span(sprintf("(%.3f)", nn.dist[1, .t]), style="color: #999999; font-style: italic;"))
        })
      do.call(tags$ol, list.items)
    }
  })

  observe({
    click <- input$mapClick
    coord <- isolate(ENV$coord)
    if (!is.null(coord) && !is.null(click)) {
      click.coord <- cbind(click$x, click$y)
      scale <- diff(range(coord[, 1])) # coordinate range on x-axis
      distances <- dist.matrix(click.coord, coord, method="euclidean")
      idx <- which.min(distances)
      if (distances[idx] < scale/25) {
        target <- rownames(coord)[idx]
        updateSelectizeInput(session, "target", choices=TargetInfo, selected=target, server=TRUE) # triggers re-display
      }
    }
  })

  output$debug <- renderPrint({
    cat(sprintf("n = %d neighbours for target: %s\n\n", input$NN, input$target))
    if (!is.null(ENV$coord)) {
      print(ENV$coord)
    }
  })

})
