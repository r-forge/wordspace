\name{t.dsm}
\alias{t.dsm}
\title{
  Swap the Rows and Columns of a DSM Object (wordspace)
}
\description{

  Given a distributional model \code{x}, \code{t(x)} returns a new DSM object
  representing the transposed co-occurrence and/or score matrix.  Marginal frequencies
  and other row/column information are updated accordingly.
  
}
\usage{

\method{t}{dsm}(x)

}
\arguments{

  \item{x}{
    an object of class \code{dsm}
  }

}
%% \details{
%% }
\value{
  An \code{DSM} object with rows and columns swapped.
}
\seealso{
  TODO
}
\examples{
## TODO
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
\keyword{TODO}
