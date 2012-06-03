\name{read.dsm.triplet}
\alias{read.dsm.triplet}
\title{
  Load DSM Data from Triplet Representation (wordspace)
}
\description{

  This functions loads DSM in triplet representation -- (row label, column label, score) -- from a disk file or a pipe.
  Such a triplet file normally represents a pre-scored DSM, but can also hold raw syntactic cooccurrence frequencies for an 
  unreduced set of target and feature terms (from which marginals and sample size can then be derived).

}
\usage{

  read.dsm.triplet(filename, freq = FALSE, value.first = FALSE, sep = "\t", quote = "",
                   nmax = -1, sort = FALSE, encoding = getOption("encoding"), verbose = FALSE) 

}
\arguments{

  \item{filename}{
    the name of a file containing the triplet data (see \sQuote{File Format} below for details), which may be compressed (\file{.gz}, \file{.bz2}, \file{.xz}).  If \code{filename} ends in \code{|}, it is opened as a Unix pipe for reading.
  }

  \item{freq}{
  
  }

  \item{value.first}{
  
  }

  \item{sep, quote}{
  
  }

  \item{nmax}{
  
  }

  \item{sort}{
  
  }

  \item{encoding}{
    character encoding of the input files, which will automatically be converted to R's internal representation if possible.  See \sQuote{Encoding} in \link{file} for details.
  }

  \item{verbose}{
    if \code{TRUE}, a few progress and information messages are shown
  }

}
\details{

  The \code{encoding} argument currently doesn't work (as of R 2.12.1) if the DSM data are packaged in a ZIP archive.  If you need to work with non-ASCII row/column labels in R (or want them to be displayed properly), unpack the archive first in a subdirectory and specify this directory as \code{filename}.  Keep in mind that individual files may still be compressed with \code{gzip} in order to save disk space.

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
