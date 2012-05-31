# this function reads the contents of a UCS export directory and returns the necessary information for .access.file()
.open.archive <- function (filename) {
  if (file.access(filename, 4) == 0) {
    if (file.info(filename)$isdir) {
        contents <- list.files(filename)
        return(list(contents=contents, dir=filename))
    } else {
      stop("file '", filename, "' exists, but is not a directory")
    }
  }
  else {
    stop("UCS export archive '", filename, "' does not exist")
  }
}

# access a file from UCS export archive, returning a connection object opened in "rt" mode;
# this function automatically finds compressed files with extension .gz, .bz2 or .xz
.access.file <- function (archive, member, encoding=getOption("encoding"), check.only=FALSE) {
  if (grepl("\\.(gz|bz2|xz)$", member, ignore.case=TRUE, perl=TRUE)) stop("internal error - archive members can only be accessed by uncompressed name")
  member.names <- paste(member, c("", ".gz", ".bz2", ".xz"), sep="")
  idx.found <- member.names %in% archive$contents
  n.found <- sum(idx.found)
  
  if (check.only) return(n.found == 1)

  if (n.found > 1) stop("UCS export archive '", archive$dir, "' contains multiple versions of the same component: ", paste(member.names[idx.found], collapse=", "))
  if (n.found == 0) stop("'", archive$dir, "' is not a valid UCS export archive, component '", member, "' missing")

  file(paste(archive$dir, member.names[idx.found], sep="/"), encoding=encoding, open="rt")
}

read.dsm.ucs <- function (filename, encoding=getOption("encoding")) {
  archive <- .open.archive(filename)
  
  fh <- .access.file(archive, "rows.tbl", encoding)
  rows <- read.delim(fh, as.is=TRUE, quote="")
  close(fh)
  
  fh <- .access.file(archive, "cols.tbl", encoding)
  cols <- read.delim(fh, as.is=TRUE, quote="")
  close(fh)
  
  n.rows <- nrow(rows)
  n.cols <- nrow(cols)
  rows$f <- as.double(rows$f) # adjust marginal frequencies to doubles
  cols$f <- as.double(cols$f)

  fh <- .access.file(archive, "globals.tbl", encoding)
  globals <- read.delim(fh, colClasses=c(N="double"), quote="")
  close(fh)

  if (nrow(globals) != 1) stop("format error - globals.tbl must contain exactly one row")
  if (!("N" %in% colnames(globals))) stop("format error - sample size N missing from globals.tbl")
  N <- globals$N <- as.double(globals$N)
  
  have.dense.M <- .access.file(archive, "M", check.only=TRUE)
  have.sparse.M <- .access.file(archive, "M.mtx", check.only=TRUE)
  if (have.dense.M && have.sparse.M) stop("UCS export archive '", archive$dir, "' contains both sparse (M.mtx) and dense (M) cooccurrence matrix")
  
  if (have.sparse.M) {
    fh <- .access.file(archive, "M.mtx", encoding)
    M <- as(readMM(fh), "dgCMatrix") # dgCMatrix (column-compressed) is the preferred format in the Matrix package
    close(fh)
  } else {
    fh <- .access.file(archive, "M", encoding)
    n.cells <- as.double(n.rows) * n.cols # avoid integer overflow for oversized matrix
    if (n.cells >= 2^31) stop("dense co-occurrence matrix is too large to load into R")
    tmp <- scan(fh, what=double(0), nmax=n.cells, quiet=TRUE)
    close(fh)
    if (length(tmp) != n.cells) stop("invalid data - M does not contain exactly ", n.cells, " = ", n.rows, " * ", n.cols, " cells")
    M <- matrix(tmp, nrow=n.rows, ncol=n.cols, byrow=TRUE)
    rm(tmp)
  }
  
  stopifnot(nrow(M) == n.rows && ncol(M) == n.cols)
  rownames(M) <- rows$term
  colnames(M) <- cols$term

  is.nzero <- M > 0 # if M is sparse, this returns a sparse logical matrix
  rows$nnzero <- rowSums(is.nzero)
  cols$nnzero <- colSums(is.nzero)
  
  dsm <- list(M=M, rows=rows, cols=cols, N=N, globals=globals, locked=FALSE)
  class(dsm) <- c("dsm", "list")
  return(dsm)
}
