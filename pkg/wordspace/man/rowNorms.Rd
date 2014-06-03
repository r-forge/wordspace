\name{rowNorms}
\alias{rowNorms}
\alias{colNorms}
\title{
  Compute Norms of Row and Column Vectors of a Matrix (wordspace)
}
\description{

  Efficiently compute the norms of all row or column vectors of a dense or sparse matrix.

}
\usage{

rowNorms(M, method = c("euclidean", "maximum", "manhattan", "minkowski"), p = 2)

colNorms(M, method = c("euclidean", "maximum", "manhattan", "minkowski"), p = 2)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{M}{
    a dense or sparse numeric matrix
  }

  \item{method}{
    norm to be computed (see \dQuote{Norms} below for details)
  }

  \item{p}{
    exponent of the \code{minkowski} p-norm, a numeric value in the range \eqn{1 \le p \le \infty}
  }

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\section{Norms}{

  \bold{TODO}
  
}
\value{

  A numeric vector containing one norm value for each row or column of \code{M}.

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
