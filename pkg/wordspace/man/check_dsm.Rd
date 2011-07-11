\name{check.dsm}
\alias{check.dsm}
\title{
  Validate Internal Structure of DSM Object (wordspace)
}
\description{

  %%  ~~ A concise (1-5 lines) description of what the function does. ~~

}
\usage{

 check.dsm(model, validate = FALSE, check.S = TRUE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{model}{
  }
  
  \item{validate}{carry out extended validation of internal consistency?}

  \item{check.S}{if set to \code{FALSE}, raw matrix \code{M} will be used as main model matrix even if scored matrix \code{S} is available (see below)}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

  Aborts with error message if any inconsistency is detected.
  Otherwise a list with the following elements is returned:

  \item{nrow}{number of rows (target terms) of the DSM}
  \item{ncol}{number of columns (features) of the DSM}
  \item{N}{sample size of the underlying data set}
  \item{slots}{character vector naming the components of the DSM object}
  \item{have.S}{whether a scored/weighted matrix S is available}
  \item{locked}{whether matrix combines data with inconsistent row or column marginals (in this case, association scores cannot be computed any more)}
  \item{sparse}{whether the main model matrix (either \code{M} or \code{S}) is stored in sparse representation or not (always uses \code{M} if \code{check.S=FALSE})}
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
