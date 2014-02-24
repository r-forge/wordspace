\name{normalize.rows}
\alias{normalize.rows}
\title{
  Normalize Rows of Matrix to Unit Length (wordspace)
}
\description{

  %%  ~~ A concise (1-5 lines) description of what the function does. ~~

}
\usage{

normalize.rows(M, \dots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{M}{
    a dense or sparse numeric matrix
  }

  \item{\dots}{
    optional \code{method} and \code{p} arguments are passed to \code{\link{rowNorms}}
  }

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

  A row-normalized matrix with the same dimensions as \eqn{M}.

}
\references{
%% ~put references to the literature/web site here ~
}
\seealso{
  See \code{\link{rowNorms}} for details on available norms and their parameters.
}
\examples{
## TODO
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
\keyword{TODO}
