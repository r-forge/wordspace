\name{merge.dsm}
\alias{merge.dsm}
\title{
  Merge Rows or Columns from Different DSM Objects (wordspace)
}
\description{

  %%  ~~ A concise (1-5 lines) description of what the function does. ~~

}
\usage{

  \method{merge}{dsm}(x, y, \dots, rows=TRUE, all=FALSE, term.suffix=NULL)

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

  \item{term.suffix}{
    optional character vector specifying one suffix string for each DSM, which will be appended to row (\code{rows=TRUE}) or column (\code{rows=FALSE}) labels in order to make them unique
  }

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  if \code{term.suffix} is specified, row information of returned DSM object will be extended with variables \code{orig.term} specifying the original terms and \code{orig.part} specifying the original component model (identified by the corresponding entry from \code{term.suffix}) 

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
