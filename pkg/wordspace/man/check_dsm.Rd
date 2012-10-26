\name{check.dsm}
\alias{check.dsm}
\title{
  Validate Internal Structure of DSM Object (wordspace)
}
\description{

  %%  ~~ A concise (1-5 lines) description of what the function does. ~~

}
\usage{

 check.dsm(model, validate = FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{model}{
  }
  
  \item{validate}{carry out extended validation of internal consistency?}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

  Aborts with error message if any inconsistency is detected.
  Otherwise a list with the following elements is returned:

  \item{nrow}{number of rows (target terms) of the DSM}
  \item{ncol}{number of columns (features) of the DSM}
  \item{N}{sample size of the underlying data set (may be \code{NA})}
  \item{slots}{character vector naming the components of the DSM object}
  \item{have.M}{whether a raw frequency matrix M is available}
  \item{have.S}{whether a score matrix S is available}
  \item{locked}{\code{TRUE} if matrix combines data with inconsistent row or column marginals (in this case, association scores cannot be computed any more)}
  \item{sparse}{whether the model matrix (\code{M} or \code{S}, or both) is stored in sparse format or not}
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
