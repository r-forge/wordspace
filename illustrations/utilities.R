##
##  Collection of utility functions for illustrations
##

library(wordspace)

# draw arrow with start/end shortened by configurable amount and optional rotated label
draw.arrow <- function (x1, y1, x2, y2, label=NULL, label.pos=0.5, label.offset=5, cut1=0, cut2=0, head1=FALSE, head2=TRUE, length=.15, lwd=1, ...) {
  .dx <- x2 - x1
  .dy <- y2 - y1
  .l <- sqrt(.dx^2 + .dy^2)
  .code <- head2 * 1 + head1 * 2
  arrows(x1 + cut1 * .dx/.l, y1 + cut1 * .dy/.l, x2 - cut2 * .dx/.l, y2 - cut2 * .dy/.l, code=.code, length=length, lwd=lwd, ...)
  if (!missing(label)) {
    .tx <- x1 + label.pos * .dx - label.offset * .dy/.l
    .ty <- y1 + label.pos * .dy + label.offset * .dx/.l
    .angle <- 180 * atan2(.dy, .dx) / pi
    .angle <- ((.angle + 90) %% 180) - 90
    text(.tx, .ty, labels=label, srt=.angle, ...)
  }
}

## -- TODO: general arguments for these functions --
# find nearest neighbours of specified noun (lemma), with optional visualisation
library(MASS)
neighbours <- function (M, term, n=10, plot=c("none", "pca", "mds", "sammon"), n.plot=n, aspect=1, iterate=FALSE, normalised=FALSE, edges=FALSE, edges.range=NULL, edges.lwd=6, edges.col="#AABBFF") {
  plot <- match.arg(plot)
  auto.edgerange <- missing(edges.range) # so range is automatically recomputed in the iterative plot
  if (! term %in% rownames(M)) stop("Noun '", term, "' is not in database.")
  term.vec <- M[term,,drop=FALSE] # make sure vector is interpreted as single-row matrix
  d <- cosine(M, term.vec, normalised=normalised, angles=TRUE)
  idx <- order(d)
  idx.show <- idx[2:(n+1)]  # first "neighbour" is term itself at distance 0
  if (plot != "none") {
    if (n.plot < 5) iterate <- FALSE
    iterate.n <- if (iterate) 5:n.plot else n.plot
    coords <- NULL
    for (n.plot in iterate.n) {
      idx.plot <- idx[1:(n.plot+1)] # map term vector and n.plot nearest neighbours into 2D plane
      N <- M[idx.plot, ]
      N.d <- cosine(N, normalised=normalised, angles=TRUE) # distance matrix (for MDS and similarity graph)
      if (plot == "pca") {
        P <- prcomp(N, center=TRUE, scale=FALSE)
        coords <- P$x[, 1:2]
      } else {
        if (iterate && !is.null(coords)) {
          # add dummy entry for new data point, then shake well
          coords <- jitter(rbind(coords, c(0,0)), amount=10) 
          coords <- sammon(N.d, y=coords, k=2, trace=FALSE)$points
        } else {
          if (plot == "sammon" || iterate) {
            coords <- sammon(N.d, k=2, trace=FALSE)$points
          }
          else {
            coords <- isoMDS(N.d, k=2, p=2, trace=FALSE)$points
          }
        }
      }
      labels <- rownames(M)[idx.plot]
      x.range <- extendrange(coords[, 1], f=.05)
      y.range <- extendrange(coords[, 2], f=.05)
      .asp <- diff(x.range) / diff(y.range)
      if (.asp < aspect) {
        x.range <- extendrange(x.range, f=(aspect/.asp-1)/2)
      } else if (.asp > aspect) {
        y.range <- extendrange(y.range, f=(.asp/aspect-1)/2)
      }
      if (iterate) Sys.sleep(.2)
      plot(coords, pch=20, xlim=x.range, ylim=y.range, xlab="", ylab="")
      if (edges) {
        ## draw connecting lines between points
        .midx <- t(combn(1:nrow(N), 2)) # 2d-index for upper triangle of distance matrix N.d
        .P1 <- .midx[,1] # idx of start node of edge
        .P2 <- .midx[,2] # idx of end node of edge
        .angle <- N.d[.midx] # angle between P1 and P2
        if (auto.edgerange) {
          edges.range <- c(mean(.angle) - 3 * sd(.angle), mean(.angle)) # show half of edges
          cat(sprintf("[showing edges for range %.1f ... %.1f degrees]\n", edges.range[1], edges.range[2]))
        }
        .lwd <- edges.lwd * (edges.range[2] - .angle) / diff(edges.range) # go from 0 to max. lwd in specified range
        .lwd <- pmin(.lwd, edges.lwd)
        .draw <- .lwd > 0.1 # draw only edges with minimum lwd of 0.1
        segments(coords[.P1[.draw],1], coords[.P1[.draw],2], coords[.P2[.draw],1], coords[.P2[.draw],2], lwd=.lwd[.draw], col=edges.col)
        points(coords, pch=20) # redraw points covered by edges
      }
      points(coords[1,1], coords[1,2], pch=20, cex=1.2, col="red")
      text(coords[-1, ], labels=labels[-1], pos=3, font=2)
      text(coords[1,1], coords[1,2], labels=labels[1], pos=3, font=2, cex=1.2, col="red")
    }
  }
  round(d[idx.show, ], 2) 
}

