\name{dsm.score}
\alias{dsm.score}
\encoding{utf8}
\title{
  Weighting, Scaling and Normalisation of Co-occurrence Matrix (wordspace)
}
\description{

  Compute feature scores for a term-document or term-term co-occurrence matrix, using one of several standard association measures.  Scores can optionally be rescaled with an isotonic transformation function and centered or standardized.  In addition, row vectors can be normalized to unit length wrt. a given norm. 

  This function has been optimized for efficiency and low memory overhead.
  
}
\usage{

dsm.score(model, score = "frequency",
          sparse = TRUE, negative.ok = NA,
          transform = c("none", "log", "root", "sigmoid"),
          scale = c("none", "standardize", "center", "scale"),
          normalize = FALSE, method = "euclidean", p = 2,
          matrix.only = FALSE, update.nnzero = FALSE)

}
\arguments{
  
  \item{model}{a DSM model, i.e. an object of class \code{dsm}}

  \item{score}{the association measure to be used for feature weighting (see \dQuote{Details} below)}

  \item{sparse}{if \code{TRUE} (the default), compute sparse non-negative association scores (see \dQuote{Details} below).
    Non-sparse association scores are only allowed if \code{negative.ok=TRUE}.
  }

  \item{negative.ok}{whether operations that introduce negative values into the score matrix (non-sparse association scores, standardization of columns, etc.) are allowed.
    The default (\code{negative.ok=NA}) is \code{TRUE} if the co-occurrence matrix \eqn{M} is dense, and \code{FALSE} if it is sparse.
  }
    
  \item{transform}{scale transformation to be applied to association scores (see \dQuote{Details} below)}

  \item{scale}{if not \code{"none"}, standardize columns of the scored matrix by z-transformation (\code{"standardize"}), center them without rescaling (\code{"center"}), or scale to unit RMS without centering (\code{"scale"})}

  \item{normalize}{if \code{TRUE} normalize row vectors of scored matrix to unit length, according to the norm indicated by \code{method} and \code{p}}

  \item{method, p}{norm to be used with \code{normalize=TRUE}.
    See \code{\link{rowNorms}} for admissible values and details on the corresponding norms
  }
  
  \item{matrix.only}{whether to return updated DSM model (default) or only the matrix of scores (\code{matrix.only=TRUE})}
  
  \item{update.nnzero}{if \code{TRUE} and a full DSM model is returned, update the counts of nonzero entries in rows and columns according to the matrix of scores (there may be fewer nonzero entries with sparse association scores, or more from dense association scores and/or column scaling)}

}
\details{
  \subsection{Association measures}{
    The following association measures can be used for feature scoring. Equations are given in the notation of Evert (2008). The most important symbols are \eqn{O_{11}}{E11} for the observed co-occurrence frequency, \eqn{E_{11}}{O11} for the co-occurrence frequency expected under a null hypothesis of independence, \eqn{R_1}{R1} for the marginal frequency of the target term, \eqn{C_1}{C1} for the marginal frequency of the feature term or context, and \eqn{N} for the sample size of the underlying corpus.  Evert (2008) explains in detail how these values are computed for different types of co-occurrence.
    
    \describe{
      \item{\code{frequency} (default)}{
        Co-occurrence \bold{frequency}: \deqn{
          O_{11}
        }{ O11 }
        Use this association measure to operate on raw, unweighted co-occurrence frequency data.
      }
      \item{\code{MI}}{
        \bold{(Pointwise) Mutual Information}, a log-transformed version of the ratio between observed and expected co-occurrence frequency: \deqn{
          \log_2 \frac{O_{11}}{E_{11}}
        }{ log2(O11 / E11) }
        Pointwise MI has a very strong bias towards pairs with low expected co-occurrence frequency (because of \eqn{E_{11}}{E11} in the denominator). It should only be applied if low-frequency targets and features have been removed from the DSM.
        
        The sparse version of MI (with negative scores cut off at 0) is sometimes referred to as "positive pointwise Mutual Information" (\bold{PPMI}) in the literature.
      }
      \item{\code{simple-ll}}{
        Simple \bold{log-likelihood} (Evert 2008, p. 1225): \deqn{
          \pm 2 \left( O_{11}\cdot \log \frac{O_{11}}{E_{11}} - (O_{11} - E_{11}) \right)
        }{ ± 2 * ( O11 * log(O11 / E11) - (O11 - E11) ) }
        This measure provides a good approximation to the full log-likelihood measure (Evert 2008, p. 1235), but can be computed much more efficiently. It is also very similar to the \bold{local-MI} measure used by several popular DSMs. The implementation used here computes \emph{signed} association scores, which are negative iff \eqn{O_{11} < E_{11}}{O11 < E11}.
        
        Log-likelihood has a strong bias towards high co-occurrence frequency and often produces a highly skewed distribution of scores. It may therefore be advisable to combine it with an additional \code{log} transformation.
      }
      \item{\code{t-score}}{
        The \bold{t-score} association measure, which is popular for collocation identification in computational lexicography: \deqn{
          \frac{O_{11} - E_{11}}{\sqrt{O_{11}}}
        }{ (O11 - E11) / sqrt(O11) }
        T-score is known to filter out low-frequency data effectively.
      }
      \item{\code{z-score}}{
        The \bold{z-score} association measure, based on a normal approximation to the binomial distribution of co-occurrence by chance: \deqn{
          \frac{O_{11} - E_{11}}{\sqrt{E_{11}}}
        }{ (O11 - E11) / sqrt(E11) }
        Z-score has a strong bias towards pairs with low expected co-occurrence frequency (because of \eqn{E_{11}}{E11} in the denominator). It should only be applied if low-frequency targets and features have been removed from the DSM.
      }
      \item{\code{Dice}}{
        The \bold{Dice coefficient} of association, which corresponds to the harmonic mean of the conditional probabilities \eqn{P(\mathrm{feature}|\mathrm{target})}{P(feature | target)} and \eqn{P(\mathrm{target}|\mathrm{feature})}{P(target | feature)}: \deqn{
          \frac{2 O_{11}}{R_1 + C_1}
        }{ 2 O11 / (R1 + C1) }
        Note that Dice is inherently sparse: it preserves zeroes and does not produce negative scores.
      }
    }
    
    The following additional scoring functions can be selected:
    \describe{
      \item{\code{tf.idf}}{
        The \bold{tf-idf} weighting scheme popular in Information Retrieval: \deqn{
          O_{11}\cdot \log \frac{1}{\mathit{df}}
        }{ O11 * log(1 / df) }
        where \eqn{\mathit{df}}{df} is the relative document frequency of the corresponding feature term and should be provided as a variable \code{df} in the model's column information.  Otherwise, it is approximated by the feature's nonzero count \eqn{n_p}{np} (variable \code{nnzero}) divided by the number \eqn{K} of rows in the co-occurrence matrix: \deqn{
          \mathit{df} = \frac{n_p + 1}{K + 1}
        }{ df = (np + 1) / (K + 1) }
        The discounting avoids division-by-zero errors when \eqn{n_p = 0}{np = 0}.
      } 
      \item{\code{reweight}}{
        Apply scale transformation, column scaling and/or row normalization to previously computed feature scores (from \code{model$S}).  This is the only \code{score} that can be used with a DSM that does not contain raw co-occurrence frequency data.
      }
    }
  }
  \subsection{Sparse association scores}{
    If \code{sparse=TRUE}, negative association scores are cut off at 0 in order to (i) ensure that the scored matrix is non-negative and (ii) preserve sparseness.  The implementation assumes that association scores are always \eqn{\leq 0}{≤ 0} for \eqn{O_{11} = 0}{O11 = 0} in this case and only computes scores for nonzero entries in a sparse matrix.  All built-in association measures satisfy this criterion.
    
    Other researchers sometimes refer to such sparse scores as "positive" measures, most notably positive point-wise Mutual Information (PPMI). Since \code{sparse=TRUE} is the default setting, \code{score="MI"} actually computes the PPMI measure.
  }
  \subsection{Scale transformations}{
    Association scores can be re-scaled with an isotonic transformation function that preserves sign and ranking of the scores. This is often done in order to de-skew the distribution of scores or as an approximate binarization (presence vs. absence of features).  The following built-in transformations are available:
    \describe{
      \item{\code{none} (default)}{
        A \bold{linear} transformation leaves association scores unchanged. \deqn{
          f(x) = x
        }
      }
      \item{\code{log}}{
        The \bold{logarithmic} transformation has a strong de-skewing effect.  In order to preserve sparseness and sign of association scores, a signed and discounted version has been implemented. \deqn{
          f(x) = \mathop{\mathrm{sgn}}(x) \cdot \log (|x| + 1)
        }{ f(x) = sgn(x) * log(|x| + 1) }
      }
      \item{\code{root}}{
        The \bold{signed square root} transformation has a mild de-skewing effect. \deqn{
          f(x) = \mathop{\mathrm{sgn}}(x) \cdot \sqrt{|x|}
        }{ f(x) = sgn(x) * sqrt(|x|) }
      }
      \item{\code{sigmoid}}{
        The \bold{sigmoid} transformation produces a smooth binarization where negative values saturate at \eqn{-1}, positive values saturate at \eqn{+1} and zeroes remain unchanged. \deqn{
          f(x) = \tanh x
        }{ f(x) = tanh(x) }
      }
    }
  }
}
\value{
  Either an updated DSM model of class \code{dsm} (default) or the matrix of (scaled and normalised) association scores (\code{matrix.only=TRUE}).
  
  Note that updating DSM models may require a substantial amount of temporary memory (because of the way memory management is implemented in \R).  This can be problematic when running a 32-bit build of \R or when dealing with very large DSM models, so it may be better to return only the scored matrix in such cases.
}
\references{

  More information about assocation measures and the notation for contingency tables can be found at \url{http://www.collocations.de/} and in

  Evert, Stefan (2008). Corpora and collocations. In A. Lüdeling and M. Kytö (eds.), \emph{Corpus Linguistics. An International Handbook}, chapter 58, pages 1212--1248. Mouton de Gruyter, Berlin, New York.

}
\seealso{

  \code{\link{dsm}}
  
}
\examples{
model <- DSM_TermTerm
model$M # raw co-occurrence matrix
  
model <- dsm.score(model, score="MI")
round(model$S, 3) # PPMI scores
  
model <- dsm.score(model, score="reweight", transform="sigmoid")
round(model$S, 3) # additional sigmoid transformation
  
\dontrun{
# visualization of the scale transformations implemented by dsm.score
x <- seq(-2, 4, .025)
plot(x, x, type="l", lwd=2, xaxs="i", yaxs="i", xlab="x", ylab="f(x)")
abline(h=0, lwd=0.5); abline(v=0, lwd=0.5)
lines(x, sign(x) * log(abs(x) + 1), lwd=2, col=2)
lines(x, sign(x) * sqrt(abs(x)), lwd=2, col=3)
lines(x, tanh(x), lwd=2, col=4)
legend("topleft", inset=.05, bg="white", lwd=3, col=1:4,
       legend=c("none", "log", "root", "sigmoid"))
}
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
