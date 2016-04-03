\name{DSM_Vectors}
\alias{DSM_Vectors}
\docType{data}
\title{
  Pre-Compiled DSM Vectors for Selected Words (wordspace)
}
\description{

  A matrix of 100-dimensional pre-compiled DSM vectors for selected English content words, covering all nouns needed for several basic evaluation tasks.
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

  The vocabulary of this DSM covers several basic evaluation tasks, including \code{\link{RG65}}, \code{\link{WordSim353}} and \code{\link{ESSLLI08_Nouns}}.  In addition 40 nearest neighbours each of the words \code{white_J}, \code{apple_N}, \code{kindness_N} and \code{walk_V} are included.

  Co-occurrence frequency data were extracted from a collection of Web corpora with a total size of ca. 9 billion words, using a L4/R4 surface window and 30,000 lexical words as feature terms.  They were scored with simple log-likelihood with an additional log transformation, normalized to Euclidean unit length, and projected into 1000 latent dimensions using randomized SVD (see \code{\link{rsvd}}.  For size reasons, only the first 100 SVD dimensions are included here, and vectors have been renormalized.
  
}
\examples{

nearest.neighbours(DSM_Vectors, "walk_V", 25)

}
\keyword{datasets}
