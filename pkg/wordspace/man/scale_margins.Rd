\name{scaleMargins}
\alias{scaleMargins}
\title{
  Scale Rows and/or Columns of a Matrix (wordspace)
}
\description{

  This function provides a reasonably memory-efficient way to scale the rows and/or columns
  of a dense or sparse matrix.

}
\usage{

  scaleMargins(M, rows=NULL, cols=NULL)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{M}{
    a dense or sparse matrix.  A sparse matrix must be given in \code{\link[=dgCMatrix-class]{dgCMatrix}} format.
  }

  \item{rows}{
    a numeric vector with length equal to the number of rows of M, or a single number.
    If missing or \code{NULL}, the rows of M are not rescaled.
  }

  \item{cols}{
    a numeric vector with length equal to the number of columns of M, or a single number.
    If missing or \code{NULL}, the columns of M are not rescaled.
  }

}
\details{
  
  NB: This function should be called \code{scale.margins}, but then it would be misinterpreted as an S3 method
  for the generic function \code{\link{scale}}.

}
\value{

  The rescaled dense or sparse matrix.

}
\references{
%% ~put references to the literature/web site here ~
}
\seealso{
}
\examples{
## TODO
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
\keyword{TODO}
