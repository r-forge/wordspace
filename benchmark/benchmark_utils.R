benchmark <- function (expr, name, n.ops=NULL, silent=FALSE) {
  gc(reset=TRUE)
  time.info <- system.time(expr)
  mem.info <- gc()
  elapsed <- time.info[3]
  user <- time.info[1]
  mops <- if (missing(n.ops)) NA else (n.ops / 1e6) / elapsed
  mb.tmp <- mem.info[2, 6] - mem.info[2, 2]
  res <- data.frame(time=elapsed, MOPS=round(mops, 2), CPUTime=user, MB=mb.tmp, row.names=name)
  if (!silent) print(res)
  res
}

append.list <- function (L, element) {
  append(L, list(element))
}

vector.equal <- function(x, y, name="vector comparison", tol=1e-12, verbose=TRUE) {
  if (length(x) == length(y)) {
    max.diff <- max(abs(x - y))
    if (max.diff < tol) {
      invisible(TRUE)
    } else {
      if (verbose) cat(sprintf("%s: largest difference between vectors = %g exceeds tolerance limit\n", name, max.diff))
      invisible(FALSE)
    }
  } else {
    if (verbose) cat(sprintf("%s: different vector lengths %d != %d\n", name, length(x), length(y)))
    invisible(FALSE)
  }
}

matrix.equal <- function(x, y, name="matrix comparison", tol=1e-12, ignore.sign=FALSE, verbose=TRUE) {
  if (nrow(x) == nrow(y) && ncol(x) == ncol(y)) {
    if (ignore.sign) {
      ## sign of columns is arbitrary in SVD projection
      max.diff.col.1 <- apply(abs(x - y), 2, max) # max diff in each column (same sign)
      max.diff.col.2 <- apply(abs(x + y), 2, max) # max diff in each column (opposite sign)
      max.diff <- max(pmin(max.diff.col.1, max.diff.col.2)) # pick smaller diff for each col, then take maximum
    } else {
      max.diff <- max(abs(x - y))
    }
    if (max.diff < tol) {
      invisible(TRUE)
    } else {
      if (verbose) cat(sprintf("%s: largest difference between matrices = %g exceeds tolerance limit\n", name, max.diff))
      invisible(FALSE)
    }
  } else {
    if (verbose) cat(sprintf("%s: different matrix formats %d x %d != %d x %d\n", name, nrow(x), ncol(x), nrow(y), ncol(y)))
    invisible(FALSE)
  }
}
