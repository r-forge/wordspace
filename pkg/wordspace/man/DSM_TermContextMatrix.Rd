\name{DSM_TermContextMatrix}
\alias{DSM_TermContextMatrix}
\docType{data}
\title{
  Example of a Term-Context Co-occurrence Matrix (wordspace)
}
\description{

  This matrix is a typical example of a term-context DSM co-occurrence matrix,
  derived from the English Wikipedia.

}
\usage{

DSM_TermContextMatrix

}
\format{

  A numeric matrix with 7 rows and 7 columns.

  Rows represent the target nouns \emph{cat}, \emph{dog}, \emph{animal}, \emph{time}, \emph{reason}, \emph{cause}, \emph{effect}.

  Columns specify the occurrence frequencies of these nouns in Wikipedia articles on \emph{Felidae}, \emph{Pet}, \emph{Feral}, \emph{Boat}, \emph{Philosohpy}, \emph{Kant} and \emph{Back Pain}.
  
}
%% \details{
%% }
%% \source{
%% }
\seealso{

  This matrix describes the same target nouns as the term-term matrix \code{\link{DSM_TermTermMatrix}}.

}
\examples{

# M M' = symmetric matrix of co-occurrence frequencies of nouns within articles
tcrossprod(DSM_TermContextMatrix) 

}
\keyword{datasets}
