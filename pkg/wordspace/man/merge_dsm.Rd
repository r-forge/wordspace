\name{merge.dsm}
\alias{merge.dsm}
\title{
  Merge Rows or Columns from Different DSM Models (wordspace)
}
\description{

  %%  ~~ A concise (1-5 lines) description of what the function does. ~~

}
\usage{

  \method{merge}{dsm}(x, y, \dots, rows=TRUE, all=FALSE)

}
\arguments{

  \item{x, y, \dots}{
    two or more objects of class \code{dsm} to be merged
  }

  \item{rows}{
    whether to merge rows (TRUE, default) or columns (FALSE) of the DSM matrices
  }
  
  \item{all}{
    if FALSE (default), only features shared by all DSMs are included in the merged DSM (or target terms with \code{rows=FALSE}).  If TRUE, all features are included with missing frequency / score values replaced by zero (analogously for target terms with \code{rows=FALSE}).  \bold{This option is not implemented yet.}
  }

  \item{deparse.level}{
    ignored
  }

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
