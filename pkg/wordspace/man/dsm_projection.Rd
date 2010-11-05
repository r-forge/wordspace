\name{dsm.projection}
\alias{dsm.projection}
\title{
  Reduce Dimensionality of DSM by Subspace Projection (wordspace)
}
\description{

  %%  ~~ A concise (1-5 lines) description of what the function does. ~~

}
\usage{

  dsm.projection(model, method=c("svd"), n=NA)

}
\arguments{

  \item{model}{
    an object of class \code{dsm}, which has already been scored/weighted with \code{dsm.score}
  }

  \item{method}{
    projection method to use for dimensionality reduction (see "DETAILS" below)
  }
  
  \item{n}{
    an integer specifying the number of target dimensions.  If unspecified, the maximal number of latent dimensions is used (i.e. the minimum of the number of rows and columns of the DSM matrix)
  }

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A numeric matrix with \code{n} columns (latent dimensions) and the same number of rows as the original DSM.
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
