##
##  Shiny semantic similarity GUI -- Load wordspace package and pre-compiled DSM
##

library(shiny)
library(wordspace)

argv <- commandArgs(trailingOnly=TRUE)
if (!(length(argv) %in% 1:2)) stop("Usage:  R --no-save -e 'shiny::runApp()' model.rda [n.dims]")

model.file <- argv[1]
n.dim <- if (length(argv) >= 2) as.integer(argv[2]) else Inf

cat(sprintf("Loading pre-compiled DSM from %s\n", model.file))
name <- load(model.file)
M <- get(name)
cat(sprintf("%d x %d matrix (%.1f MB)\n", nrow(M), ncol(M), object.size(M) / 2^20))

rm(list=name)
gc()
if (n.dim < ncol(M)) {
  M <- M[, 1:n.dim]
  gc()
  cat(sprintf("Reduced to %d x %d matrix (%.1f MB)\n", nrow(M), ncol(M), object.size(M) / 2^20))
}
