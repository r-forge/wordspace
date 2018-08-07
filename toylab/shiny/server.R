##
##  Shiny semantic similarity GUI -- Server implementation
##

max.cand <- 100 # max. number of completion candidates to show in UI

shinyServer(function(input, output, session) {
  ENV <- reactiveValues(
    target="",
    model=NULL,
    terms=NULL,
    coord=NULL
  )
  
  observeEvent(input$candidates, {
    if (length(input$candidates) > 0) {
      candidate <- input$candidates[1]
      updateTextInput(session, "target", value=candidate) # triggers display
    }    
  })
  
  observeEvent(c(input$target, input$candidates), {
    target <- input$target
    terms.sorted <- ENV$terms
    if (grepl("%+$", target, perl=TRUE)) {
      ## explicit prefix search
      target <- sub("%+$", "", target, perl=TRUE)
      target.complete <- FALSE
    } else {
      target.complete <- target %in% terms.sorted
    }

    if (target.complete) {
      ENV$target <- target
      updateSelectInput(session, "candidates", selected=character(0))
    } else {
      if (target == "") {
        candidates <- character(0)
      } else {
        idx <- substr(terms.sorted, 1, nchar(target)) == target
        candidates <- terms.sorted[idx]
        if (length(candidates) > max.cand) candidates <- candidates[1:max.cand]
      }
      updateSelectInput(session, "candidates", choices=candidates, selected=character(0))
      if (length(candidates) == 0) ENV$target <- target # not in DSM
    } 
  })

  observeEvent(input$model, {
    name <- input$model
    if (! (name %in% names(Models))) {
      showModal(modalDialog(
        title="Loading model ...", footer=NULL,
        strong(name), sprintf("– %.1f MiB", file.info(Files[name])$size / (2^20)),
        br(), em(Files[name])
      ))
      ENV$model <- fetch.model(name)
      removeModal()
    } else {
      ENV$model <- Models[[ name ]]
    }
    ENV$terms <- Terms[[ name ]]
    updateSelectInput(session, "candidates", choices=character(0), selected=character(0))
  })

  NN <- reactive({
    target <- ENV$target
    M <- ENV$model
    ok <- target %in% ENV$terms
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
      text(1, 1, sprintf("%s not found in DSM", isolate(ENV$target)))
    } else {
      ENV$coord <- plot(nn.dist, expand=.1, cex=1.2, method=if (input$isoMDS) "isomds" else "sammon", show.edges=input$edges, edges.lwd=input$edgeWidth)
    }
  })

  output$NNlist <- renderUI({
    nn.dist <- NN()
    if (is.null(nn.dist)) {
      span(tags$em(isolate(ENV$target)), "not found in DSM", style="color: red;")
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
        updateTextInput(session, "target", value=target) # triggers re-display
        updateSelectInput(session, "candidates", choices=character(0), selected=character(0))
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
