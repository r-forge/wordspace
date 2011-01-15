# this function checks whether the DSM data are stored in a directory or ZIP archive
# and returns the necessary information for .access.file()
.open.archive <- function (filename) {
  if (file.access(filename, 4) == 0) {
    if (grepl("\\.zip$", filename, ignore.case=TRUE, perl=TRUE)) {
      contents <- unzip(filename, list=TRUE)
      return(list(contents=contents$Name, filename=filename, type="zip"))
    } else 
    if (file.info(filename)$isdir) {
        contents <- list.files(filename)
        return(list(contents=contents, filename=filename, type="dir"))
    } else {
      stop("file '", filename, "' exists, but is not a directory")
    }
  }
  else {
    filename.zip <- paste(filename, "zip", sep=".")
    if (file.access(filename.zip, 4) == 0) {
      return(.open.archive(filename.zip))
    } else {
      stop("neither '", filename, "' nor '", filename.zip, "' exists")
    }
  }
}

# access a file from archive or directory, returning a connection object opened in "rt" mode;
# this function automatically finds and decompresses a gzipped version with extension .gz
.access.file <- function (archive, member, encoding=getOption("encoding")) {
  if (grepl("\\.gz$", member, ignore.case=TRUE, perl=TRUE)) stop("internal error - archive members can only be access by uncompressed name")
  member.gz <- paste(member, "gz", sep=".")
  have.plain <- member %in% archive$contents
  have.gz <- member.gz %in% archive$contents
  if (have.plain && have.gz) warning("DSM archive '", archive$filename, "' contains both '", member, "' and '", member.gz, "' -- using uncompressed version")
  if (!(have.plain || have.gz)) stop("'", archive$filename, "' is not a valid DSM archive, entry '", member, "' missing")

  fh <- NULL
  if (archive$type == "zip") {
    if (have.plain) {
      fh <- unz(archive$filename, member, encoding=encoding)
    } else {
      stop("extraction of .gz files from .zip archive not yet implemented")
    }
  } else {
    if (have.plain) {
      fh <- file(paste(archive$filename, member, sep="/"), encoding=encoding)
    } else {
      fh <- gzfile(paste(archive$filename, member.gz, sep="/"), encoding=encoding)
    }
  }

  open(fh, "rt") # open= argument to unz() opens in strange mode that doesn't work with read.table()
  fh
}

read.dsm <- function (filename, encoding=getOption("encoding")) {
  archive <- .open.archive(filename)
  
  fh <- .access.file(archive, "rows.tbl", encoding)
  rows <- read.delim(fh, as.is=TRUE, quote="")
  close(fh)
  
  fh <- .access.file(archive, "cols.tbl", encoding)
  cols <- read.delim(fh, as.is=TRUE, quote="")
  close(fh)
  
  n.rows <- nrow(rows)
  n.cols <- nrow(cols)
  n.cells <- n.rows * n.cols

  if (any(c("globals.tbl", "globals.tbl.gz") %in% archive$contents)) {
    fh <- .access.file(archive, "globals.tbl")
    globals <- read.delim(fh, colClasses=c(N="double"), quote="")
    close(fh)
    if (nrow(globals) != 1) stop("format error - globals.tbl must contain exactly one row")
    rows$f <- as.double(rows$f) # adjust marginal frequencies to doubles
    cols$f <- as.double(cols$f)
  }
  else {
    # compatibility mode: read old-style format (deprecated, will be removed in v1.0)
    if (! all(c("R1","R2") %in% colnames(rows))) stop("format error - old-style DSM archive without row marginals R1, R2")
    if (! all(c("C1","C2") %in% colnames(cols))) stop("format error - old-style DSM archive without column marginals C1, C2")
    rows$f <- as.double(rows$R1)
    cols$f <- as.double(cols$C1)
    
    N.vec <- c(rows$f + rows$R2, cols$f + cols$C2) # determine sample size from row/column marginals
    if (diff(range(N.vec)) > 1e-12) warning(sprintf("old-style marginals with inconsistent sample sizes in range %d .. %d (using maximum)", min(N.vec), max(N.vec)))
    globals <- data.frame(N=max(N.vec))
    
    rows[["R1"]] <- NULL; rows[["R2"]] <- NULL
    cols[["C1"]] <- NULL; cols[["C2"]] <- NULL
  }

  if (!("N" %in% colnames(globals))) stop("format error - sample size N missing from globals.tbl")
  N <- globals$N <- as.double(globals$N)
  
  have.dense.M <- any(c("M", "M.gz") %in% archive$contents)
  have.sparse.M <- any(c("M.mtx", "M.mtx.gz") %in% archive$contents)
  if (have.dense.M && have.sparse.M) warning("DSM archive '", archive$filename, "' contains both sparse (M.mtx) and dense (M) cooccurrence matrix -- using M.mtx")
  if (have.sparse.M) {
    fh <- .access.file(archive, "M.mtx", encoding)
    M <- as(readMM(fh), "dgCMatrix") # dgCMatrix (column-compressed) is the preferred format in the Matrix package
    close(fh)
  } else {
    fh <- .access.file(archive, "M", encoding)
    tmp <- scan(fh, what=double(0), nmax=n.cells, quiet=TRUE)
    close(fh) # not sure when we need to / are allowed to close connections
    if (length(tmp) != n.cells) stop("invalid data - M does not contain exactly ", n.cells, " = ", n.rows, " * ", n.cols, " cells")
    M <- matrix(tmp, nrow=n.rows, ncol=n.cols, byrow=TRUE)
    rm(tmp)
  }
  
  stopifnot(nrow(M) == n.rows && ncol(M) == n.cols)
  rownames(M) <- rows$term
  colnames(M) <- cols$term
  dsm <- list(M=M, rows=rows, cols=cols, N=N, globals=globals, locked=FALSE)
  class(dsm) <- c("dsm", "list")
  return(dsm)
}
