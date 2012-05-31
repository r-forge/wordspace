\name{read.dsm.ucs}
\alias{read.dsm.ucs}
\title{
  Load Raw DSM Data from Disk Files in UCS Export Format (wordspace)
}
\description{

  This function loads raw DSM data -- a cooccurrence frequency matrix and tables of marginal frequencies -- in \bold{UCS} export format. The data are read from a directory containing several text files with predefined names, which can optionally be compressed (see \sQuote{File Format} below for details).

}
\usage{

  read.dsm.ucs(filename, encoding=getOption("encoding"))

}
\arguments{

  \item{filename}{
    the name of a directory containing files with the raw DSM data.
  }

  \item{encoding}{
    character encoding of the input files, which will automatically be converted to R's internal representation if possible.  See \sQuote{Encoding} in \code{\link{file}} for details.
  }

}
\details{

\bold{TODO}

}
\value{

  An object of class \code{dsm}, a list with the following components:

  \bold{NB:} tables \code{rows} and \code{cols} contain variables \code{term}, \code{f}, and \code{nnzero} with the number of nonzero (positive) entries in the corresponding row or column; any additional variables given in the disk files are also included with their original names
  
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  %%  \item{comp1 }{Description of 'comp1'}
  %%  \item{comp2 }{Description of 'comp2'}
  %% ...

}
\section{File Format}{

TODO

a directory or ZIP archive containing the following files with the specified names

files:
\itemize{
  \item \file{M} \strong{or} \file{M.mtx}

  cooccurrence matrix (dense, plain text) or sparse matrix (MatrixMarket format)
  
  \item \file{rows.tbl}
  
  row information (labels \code{term}, marginal frequencies \code{f})
  
  \item \file{cols.tbl}
  
  column information (labels \code{term}, marginal frequencies \code{f})

  \item \file{globals.tbl}

  table with single row containing global variables; must include variable \code{N} specifying sample size
}

each individual file may be compressed and may carry an additional filename extension \code{.gz}, \code{.bz2} or \code{.xz}; \code{read.dsm.ucs} will automatically decompress the files when loading them

}
\note{
%%  ~~further notes~~
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## TODO
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
\keyword{file}
