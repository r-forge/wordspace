\name{read.dsm.triplet}
\alias{read.dsm.triplet}
\title{
  Load DSM Data from Triplet Representation (wordspace)
}
\description{

  This functions loads a sparse distributional semantic model in triplet representation -- (row label, column label, score) -- from a disk file or a pipe.  Such a triplet file normally represents a pre-scored DSM, but can also hold raw syntactic cooccurrence frequencies for an unreduced set of target and feature terms (from which marginals and sample size can then be derived).

}
\usage{

read.dsm.triplet(filename, freq = FALSE, value.first = FALSE, tokens = FALSE, sep = "\t", quote = "",
                 nmax = -1, sort = FALSE, encoding = getOption("encoding"), verbose = FALSE) 

}
\arguments{

  \item{filename}{
    the name of a file containing the triplet data (see \sQuote{File Format} below for details), which may be compressed (\file{.gz}, \file{.bz2}, \file{.xz}).  If \code{filename} ends in \code{|}, it is opened as a Unix pipe for reading.
  }

  \item{freq}{
    whether values are raw co-occurrence frequencies (\code{TRUE}) or pre-computed scores (\code{FALSE})
  }

  \item{value.first}{
    if \code{TRUE}, triplets are given as (score, row label, column label) instead of the default (row label, column label, score)
  }

  \item{tokens}{
    if \code{TRUE}, the input file contains pair \emph{tokens}, i.e. row and column labels without score/frequency values.
    Co-occurrence frequencies will automatically be calculated, but this input format should only be used
    for small samples up to a few millon tokens.
  }

  \item{sep, quote}{
    specify field separator and the types of quotes used by the disk file (see the \code{\link{scan}} documentation for details).  By default, a TAB-delimited file without quotes is assumed.
  }

  \item{nmax}{
    if the number of entries (= text lines) in the disk file is known, it can be specified here in order to make loading faster and more memory-efficient.  Caution: If \code{nmax} is smaller than the number of lines in the disk file, the extra lines will silently be discarded.
  }

  \item{sort}{
    if \code{TRUE}, the rows and columns of the co-occurrence matrix will be sorted alphabetically according to their labels (i.e. the target and feature terms); otherwise they are listed as encountered in the triplet representation
  }

  \item{encoding}{
    character encoding of the input files, which will automatically be converted to \R's internal representation if possible.  See \sQuote{Encoding} in \code{\link{file}} for details.
  }

  \item{verbose}{
    if \code{TRUE}, a few progress and information messages are shown
  }

}
\details{

  \bold{TODO}
  
}
\value{

  An object of class \code{\link{dsm}} containing a sparse DSM.

}
\section{File Format}{

  The input must be a plain-text table with two or three TAB-delimited columns and no header.

  For \code{tokens=TRUE}, each line represents a single pair token with columns
  \enumerate{
    \item target term
    \item feature term / document
  }

  For \code{tokens=FALSE}, each line represents a pair type (i.e. a unique cell of the co-occurrence matrix) with columns:
  \enumerate{
    \item target term
    \item feature term / document
    \item score \emph{or} co-occurrence frequency
  }
  If \code{value.first=TRUE}, the score entry is expected in the first column:
  \enumerate{
    \item score \emph{or} co-occurrence frequency
    \item target term
    \item feature term / document
  }
  Note that the input may contain multiple entries for the same cell, whose values will automatically be added up.
}
\seealso{

  \code{\link{dsm}}, \code{\link{read.dsm.ucs}}
  
}
%% \examples{
%% }
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
\keyword{file}
