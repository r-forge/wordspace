wordspace.openmp <- function (threads=NULL) {
  if (!is.null(threads)) {
    if (!(is.numeric(threads) && length(threads) == 1 && threads >= 1)) stop("argument threads= must be a single integer >= 1")
    .C("C_set_openmp_threads", n=as.integer(threads), DUP=FALSE, NAOK=FALSE)
  }
  omp <- .C("C_get_openmp_threads", use=integer(1), max=integer(1), DUP=FALSE)
  res <- data.frame(available=(omp$max > 0), max=omp$max, threads=omp$use, row.names="OpenMP")
  if (is.null(threads)) res else invisible(res)
}
