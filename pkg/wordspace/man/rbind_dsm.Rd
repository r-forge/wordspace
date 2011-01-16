\name{rbind.dsm}
\alias{rbind.dsm}
\alias{cbind.dsm}
\title{
  Combine DSM Objects by Rows and Columns (wordspace)
}
\description{

  %%  ~~ A concise (1-5 lines) description of what the function does. ~~

}
\usage{

  \method{rbind}{dsm}(\dots, term.suffix=NULL, deparse.level = 1)
  \method{cbind}{dsm}(\dots, term.suffix=NULL, deparse.level = 1)

}
\arguments{

  \item{\dots}{
    one or more objects of class \code{dsm}, which must have the same feature dimensions (\code{rbind}) or target terms (\code{cbind})
  }

  \item{term.suffix}{
    optional character vector specifying one suffix string for each DSM, which will be appended to row (\code{rbind}) or column (\code{cbind}) labels in order to make them unique
  }
  
  \item{deparse.level}{
    ignored
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
