##
##  Shiny semantic similarity GUI -- Load wordspace package and pre-compiled DSM
##

library(shiny)
library(wordspace)

argv <- commandArgs(trailingOnly=TRUE)
if (length(argv) < 1) stop("Usage:  R --no-save -e 'shiny::runApp()' --args [name=]model.rda [[name2=]model2.rda ...]")

name.file <- strsplit(argv, "=")
idx <- !(sapply(name.file, length) %in% c(1, 2))
if (any(idx)) stop("Invalid argument(s): ", paste(argv[idx], collapse=", "))

Files <- sapply(name.file, function (x) x[length(x)])
Names <- sapply(name.file, function (x) x[1])
names(Files) <- Names

idx <- file.access(Files, mode=4) != 0
if (any(idx)) stop("Can't access file(s): ", paste(Files[idx], collapse=", "))

Models <- list()
Terms <- list()

fetch.model <- function (name) {
  stopifnot(name %in% Names)
  if (! (name %in% names(Models))) {
    cat(sprintf("Loading pre-compiled DSM from %s\n", Files[name]))
    varname <- load(Files[name])
    Models[[ name ]] <<- x <- get(varname)
    Terms[[ name ]] <<- sort(rownames(x))
    cat(sprintf("%d x %d matrix (%.1f MB)\n", nrow(x), ncol(x), object.size(x) / 2^20))
  }
  invisible(Models[[ name ]])
}

## fetch.model(Names[1]) # pre-load first (or only) DSM
