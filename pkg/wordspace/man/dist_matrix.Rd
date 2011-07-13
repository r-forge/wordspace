\name{dist.matrix}
\alias{dist.matrix}
\title{
  Distances/Similarities between Row or Column Vectors (wordspace)
}
\description{

  Compute a symmetric matrix of distances (or similarities) between the rows or columns of a matrix;
  or compute cross-distances between the rows or columns of two different matrices.
  This implementation is faster than \code{\link{dist}} and can operate on sparse matrices (class \code{dgCMatrix}).

}
\usage{

dist.matrix(M, M2 = NULL, method = c("cosine", "euclidean", "maximum", "manhattan", "minkowski"), p = 2, 
            normalized = FALSE, byrow = TRUE, convert = TRUE, as.dist = FALSE)

}
\arguments{

  \item{M}{
    a dense or sparse matrix representing a scored DSM
  }
  
  \item{M2}{
    an optional dense or sparse matrix representing a second scored DSM.
    If present, cross-distances between the rows (or columns) of \code{M} and those of \code{M2} will be computed.
  }

  \item{method}{
    distance or similarity measure to be used (see "Distance Measures" below for details)
  }
  
  \item{p}{
    exponent of the \code{minkowski} p-metric, a numeric value in the range \eqn{1\le p\le \infty}
  }
  
  \item{normalized}{
    if \code{TRUE}, assume that the row (or column) vectors of \code{M} and \code{M2} have been appropriately normalised (depending on the selected distance measure) in order to speed up calculations.
    This option is often used with the \code{cosine} measure, for which vectors must be normalized wrt. the Euclidean norm.
  }

  \item{byrow}{
    whether to calculate distances between row vectors (default) or between column vectors (\code{byrow=FALSE})
  }
  
  \item{convert}{
    if \code{TRUE}, similarity measures are automatically converted to distances in an appropriate way (see "Distance Measures" below for details).  Note that this is the default setting and \code{convert=FALSE} has to be specified in order to obtain a similarity matrix.
  }
  
  \item{as.dist}{
    convert the full symmetric distance matrix to a compact object of class \code{dist}.  This option cannot be used if cross-distances are calculated (i.e. if \code{M2} is specified) or if a similarity measure has been selected (with option \code{convert=FALSE}).
  }

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  By default, a matrix of class \code{dist.matrix}, specifying distances or similarities between term vectors.
  A similarity matrix is marked by an additional attribute \code{similarity} with value \code{TRUE}.
  
  If \code{as.dist=TRUE}, the matrix is compacted to an object of class \code{dist}.

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\section{Distance Measures}{
  metric vs. distance measure vs. similarity measure
  
  explain for each similarity measure how the conversion to a distance is accomplished
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
