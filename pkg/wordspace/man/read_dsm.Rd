\name{read.dsm}
\alias{read.dsm}
\title{
  Load raw DSM data from disk file
}
\description{
  This functions loads raw DSM data - a cooccurrence frequency matrix and tables of marginal frequencies - from a disk file in \emph{\pkg{wordspace}} text format.  The input is either a ZIP archive or a directory containing text files with predefined names (see \sQuote{File Format} below for details).
}
\usage{
read.dsm(filename, encoding=getOption("encoding"))
}
\arguments{
  \item{filename}{the name of a directory or ZIP archive containing the raw DSM data.  If the specified file does not exist, looks for a corresponding ZIP archive with the additional extension \code{.zip}.}
  \item{encoding}{character encoding of the input files, which will automatically be converted to R's internal representation if possible.  See \sQuote{Encoding} in \link{file} for details.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class \code{dsm}, a list with the following components:

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\section{File Format}{
TODO

either directory or ZIP archive containing the following files with the specified names

note that in the ZIP archive case, there must not be any prefix before the filenames; most GUI tools for creating ZIP archives will not handle this properly

files:
\itemize{
  \item \file{M} \strong{or} \file{M.mtx}

  cooccurrence matrix (dense, plain text) or sparse matrix (MatrixMarket format)
  
  \item \file{rows.tbl}
  
  row information (labels, marginal frequencies R1, R2)
  
  \item \file{cols.tbl}
  
  column information (labels, marginal frequencies C1, C2)
}

in the case of a directory, each contained file may be individually GZip-compressed with the extension \code{.gz}; \code{read.dsm} will automatically decompress the files when loading them

it makes little sense to compress the files within a ZIP archive and this is currently not supported by \code{read.dsm}

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
