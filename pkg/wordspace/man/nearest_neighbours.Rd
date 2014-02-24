\name{nearest.neighbours}
\alias{nearest.neighbours}
\alias{nearest.neighbors}
\title{
  Find Nearest Neighbours in DSM Space (wordspace)
}
\description{

  Find the nearest neighbours of a term vector in a DSM, given either as a scored cooccurrence matrix or a pre-computed distance matrix.  The target term is selected by name, so the cooccurrence or distance matrix must be labelled appropriately.

}
\usage{

nearest.neighbours(M, term, n = 10, drop = TRUE, skip.missing = FALSE,
                   min.distance = 1e-6, byrow = TRUE, \dots,
                   batchsize=50e6, verbose=FALSE)

}
\arguments{

  \item{M}{
    either a dense or sparse matrix representing a scored DSM (or an object of class \code{dsm}), or a pre-computed distance matrix returned by \code{dist.matrix} (as an object of class \code{dist.matrix}).  Note that the compact representation returned by the \code{dist} function (class \code{dist}) is not accepted.
  }

  \item{term}{
    a character vector specifying one or more target terms for which nearest neighbours will be found
  }

  \item{n}{
    an integer giving the number of nearest neighbours to be returned for each target term
  }

  \item{drop}{
    if \code{TRUE}, the return value is simplified to a vector if it contains nearest neighbours for exactly one target term (default).  Set \code{drop=FALSE} to ensure that \code{nearest.neighbours} always returns a list.
  }

  \item{skip.missing}{
    if \code{TRUE}, silently ignores target terms not found in the DSM or distance matrix.  By default (\code{skip.missing=FALSE}) an error is raised in this case.
  }
  
  \item{min.distance}{
    in most cases, the nearest neighbour of a target term will be the term itself at a distance of 0.  By default, this entry is automatically removed from the list of neighbours, using a threshold of \code{min.distance} to allow for rounding errors in the distance computation.  Set \code{min.distance=0} to include the target term among the nearest neighbours.  Note that this feature does \emph{not} apply to similarity measures.
  }

  \item{byrow}{
    wether target terms are looked up in rows (default) or columns (\code{byrow=FALSE}) of \code{M}
  }

  \item{\dots}{
    additional arguments are passed to \code{dist.matrix} if \code{M} is a scored DSM matrix.  See the manpage of \code{\link{dist.matrix}} for details on available parameters and settings.
  }
    
  \item{batchsize}{if \code{term} is a long list of lookup terms, it will automatically be processed in batches.  The number of terms per batch is chosen in such a way that approximately \code{batchsize} intermediate similarity values have to be computed and stored at a time.}
  
  \item{verbose}{if \code{TRUE}, display some progress messages indicating how data are split into batches}
   
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

  A list with one entry for each target \code{term} found in \code{M}, giving the nearest neighbours as a numeric vector of distances or similarities labelled with the corresponding terms and ordered by distance.

  If \code{drop=TRUE}, a list containing only a single target term will be simplified to a plain vector.

}
\references{
%% ~put references to the literature/web site here ~
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## TODO
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
\keyword{TODO}
