\name{DSM_Vectors}
\alias{DSM_Vectors}
\docType{data}
\title{
  Pre-Compiled DSM Vectors for Selected Words (wordspace)
}
\description{

  A matrix of 100-dimensional pre-compiled DSM vectors for selected content words, covering all nouns needed for several basic evaluation tasks.
  Targets are given as disambiguated lemmas in the form \code{<headword>_<pos>}, e.g. \code{walk_V} and \code{walk_N}.

}
\usage{
  
DSM_Vectors

}
\format{
  A numeric matrix with 1331 rows and 100 columns.
  
  Row labels are disambiguated lemmas of the form \code{<headword>_<pos>}, where the part-of-speech code is one of
  \code{N} (noun), \code{V} (verb), \code{J} (adjective) or \code{R} (adverb).
  
}
\details{

  \bold{TODO}: describe model compilation
  
}
\examples{

nearest.neighbours(DSM_Vectors, "walk_V", 25)

}
\keyword{datasets}
