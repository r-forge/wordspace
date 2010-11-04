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

# access a file from archive or directory, returning a connection object;
# this function automatically finds and decompresses a gzipped version with extension .gz
.access.file <- function (archive, member, encoding=getOption("encoding")) {
  if (grepl("\\.gz$", member, ignore.case=TRUE, perl=TRUE)) stop("internal error - archive members can only be access by uncompressed name")
  member.gz <- paste(member, "gz", sep=".")
  have.plain <- member %in% archive$contents
  have.gz <- member.gz %in% archive$contents
  if (have.plain && have.gz) warning("DSM archive '", archive$filename, "' contains both '", member, "' and '", member.gz, "' -- using uncompressed version")
  if (!(have.plain || have.gz)) stop("'", archive$filename, "' is not a valid DSM archive, entry '", member, "' missing")
  
  if (archive$type == "zip") {
    if (have.plain) {
      return(unz(archive$filename, member, encoding=encoding))
    } else {
      stop("extraction of .gz files from .zip archive not yet implemented")
    }
  } else {
    if (have.plain) {
      return(file(paste(archive$filename, member, sep="/"), encoding=encoding))
    } else {
      return(gzfile(paste(archive$filename, member.gz, sep="/"), encoding=encoding))
    }
  }
}

read.dsm <- function (filename, encoding=getOption("encoding")) {
  archive <- .open.archive(filename)
  
  fh <- .access.file(archive, "rows.tbl", encoding)
  rows <- read.delim(fh, colClasses=c("character", "double", "double"), quote="")

  fh <- .access.file(archive, "cols.tbl", encoding)
  cols <- read.delim(fh, colClasses=c("character", "double", "double"), quote="")

  n.rows <- nrow(rows)
  n.cols <- nrow(cols)
  n.cells <- n.rows * n.cols

  have.dense.M <- any(c("M", "M.gz") %in% archive$contents)
  have.sparse.M <- any(c("M.mtx", "M.mtx.gz") %in% archive$contents)
  if (have.dense.M && have.sparse.M) warning("DSM archive '", archive$filename, "' contains both sparse (M.mtx) and dense (M) cooccurrence matrix -- using M.mtx")
  if (have.sparse.M) {
    fh <- .access.file(archive, "M.mtx", encoding)
    M <- readMM(fh)
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
  dsm <- list(M=M, rows=rows, cols=cols, locked=FALSE)
  class(dsm) <- c("dsm", "list")
  return(dsm)
}