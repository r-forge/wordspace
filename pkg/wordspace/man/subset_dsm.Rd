\name{subset.dsm}
\alias{subset.dsm}
\title{
  Subsetting Distributional Semantic Models (wordspace)
}
\description{

  Filter the rows and/or columns of a DSM object according to user-specified conditions.

}
\usage{

\method{subset}{dsm}(x, subset, select, recursive = FALSE, drop.zeroes = FALSE, matrix.only = FALSE, envir=parent.frame(), \dots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{x}{
    an object of class \code{dsm}
  }

  \item{subset}{
    Boolean expression or index vector selecting a subset of the rows; expression can use variables \code{term} and \code{f} to access target terms and their marginal frequencies, \code{nnzero} for the number of nonzero elements in each row, further optional variables from the row information table, as well as global variables such as the sample size \code{N}
  }

  \item{select}{
    Boolean expression or index vector selecting a subset of the columns; expression can use variables \code{term} and \code{f} to access feature terms and their marginal frequencies, \code{nnzero} for the number of nonzero elements in each column, further optional variables from the column information table, as well as global variables such as the sample size \code{N}
  }

  \item{recursive}{
    if \code{TRUE} and both \code{subset} and \code{select} conditions are specified, the \code{subset} is applied repeatedly until the DSM no longer changes
  }
  
  \item{drop.zeroes}{
    if \code{TRUE}, all rows and columns without nonzero entries after subsetting are removed from the model
    (nonzero counts are based on score matrix \eqn{S} if available, raw cooccurrence frequencies \eqn{M} otherwise)
  }

  \item{matrix.only}{
    if \code{TRUE}, return only the selected subset of the score matrix \eqn{S} (if available) or frequency matrix \eqn{M}, not a full DSM object.  This may conserve a substantial amount of memory when processing very large DSMs.
  }

  \item{envir}{
    environment in which the \code{subset} and \code{select} conditions are evaluated.  Defaults to the context of the function call, so all variables visible there can be used in the expressions.
  }
  
  \item{\dots}{
    any further arguments are silently ignored
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

  An object of class \code{dsm} containing the specified subset of the model \code{x}.

  If necessary, counts of nonzero elements for each row and/or column are updated automatically.

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
