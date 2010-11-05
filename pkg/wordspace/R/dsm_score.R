dsm.score <- function (model,
                       score=c("frequency", "t-score", "z-score", "Dice", "MI", "reweight"), sparse=FALSE,
                       transform=c("none", "log", "root", "sigmoid"),
                       scale=c("none", "standardize"),
                       normalize=FALSE, method="euclidean", p=2) {
  score <- match.arg(score)
  transform <- match.arg(transform)
  scale <- match.arg(scale)
  model.info <- check.dsm(model)
  if (model.info$locked && !(score %in% c("frequency", "reweight"))) stop("marginal frequencies have been adjusted, cannot recompute association scores")
  if (score == "reweight" && !model.info$have.S) stop("cannot use score='reweight': association scores have not been computed yet")
  if (scale == "standardize" && sparse) warning("standardization of matrix columns destroys sparse representation: are you sure you want to do this?")
  if (normalize && !sparse) warning("normalization of matrix rows is only sensible for a sparse non-negative representation")
  
  O <- model$M # observed frequencies
  R1 <- model$rows$R1
  R2 <- model$rows$R2
  C1 <- model$cols$C1
  C2 <- model$cols$C2
  
  need.exp <- score %in% c("t-score", "z-score", "MI")
  if (need.exp) {
    N <- max(R1 + R2) # if there should be differences, largest value is most consistent
    if (N != max(C1 + C2)) stop("row and column marginals lead to inconsistent sample sizes N")
    E <- outer(R1, C1) / N
  }

  if (sparse) {
    scores <- matrix(0, nrow(O), ncol(O))
    nonzero <- O > 0
  }

  if (score == "t-score") {
    if (sparse) {
      idx <- nonzero & O > E
      scores[idx] <- (O[idx] - E[idx]) / sqrt(O[idx])
    } else {
      scores <- (O - E) / sqrt(O + 1) # "discounted" t-score
    }
  } else if (score == "z-score") {
    if (sparse) {
      idx <- nonzero & O > E
      scores[idx] <- (O[idx] - E[idx]) / sqrt(E[idx])
    } else {
      scores <- (O - E) / sqrt(E) # z-score (E=0 should never happen)
    }
  } else if (score == "Dice") {
    scores <- 2 * O / outer(R1, C1, "+") # Dice coefficient
  } else if (score == "MI") {
    if (sparse) {
      ids <- nonzero & O > E
      scores[idx] <- log2(O[idx] / E[idx])
    } else {
      scores <- log2(O / E) # standard pointwise MI
    }
  } else if (score == "reweight") {
    scores <- model$S # reweight existing scores (with transformation, scaling and/or normalisation)
  } else {
    scores <- O # "frequency" measure = observed frequency
  }

  if (transform == "log") {
    scores <- sign(scores) * log(abs(scores) + 1) # "signed log" transformation
  } else if (transform == "root") {
    scores <- sign(scores) * sqrt(abs(scores)) # signed square root
  } else if (transform == "sigmoid") {
    stop("not yet implemented")
  } else {
    # no transformation
  }

  rownames(scores) <- rownames(O) # make sure that row and column names are preserved
  colnames(scores) <- colnames(O)
  model$S <- scores
  return(model)
}

