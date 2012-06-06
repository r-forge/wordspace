\name{subset.dsm}
\alias{subset.dsm}
\title{
  Subsetting Distributional Semantic Models (wordspace)
}
\description{

  %%  ~~ A concise (1-5 lines) description of what the function does. ~~

}
\usage{

 \method{subset}{dsm}(x, subset, select, drop.zeroes=FALSE, matrix.only=FALSE, \dots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{x}{
  }

  \item{subset}{
    Boolean expression or index vector selecting a subset of the rows; expression can use variables \code{term} and \code{f} to access target (row) terms and their marginal frequencies, as well as optional other variables defined for the rows and global variables such as the sample size \code{N}
  }

  \item{select}{
    Boolean expression or index vector selecting a subset of the columns; expression can use variables \code{term} and \code{f} to access feature (column) terms (or context labels) and their marginal frequencies, as well as optional other variables defined for the columns and global variables such as the sample size \code{N}
  }

  \item{drop.zeroes}{
    if \code{TRUE}, all rows and columns without nonzero entries after subsetting are removed from the model
    (nonzero counts are based on score matrix S if available, raw cooccurrence frequencies M otherwise)
  }

  \item{matrix.only}{
    if \code{TRUE}, return only the selected subset of the score matrix S (if available) or frequency matrix M, not a full DSM object.
    This may conserve a substantial amount of memory when processing very large DSMs.
  }

  \item{\dots}{
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

  NB: if necessary, counts of nonzero elements in each row and/or column are updated automatically

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
