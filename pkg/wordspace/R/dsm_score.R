dsm.score <- function (model, score=c("frequency", "t-score","z-score","Dice"), sparse=FALSE, transform=c("none", "log", "root", "sigmoid")) {
  score <- match.arg(score)
  transform <- match.arg(transform)

  O <- model$M # observed frequencies
  stopifnot(nrow(O) == nrow(model$rows))
  stopifnot(ncol(O) == nrow(model$cols))
  R1 <- model$rows$R1
  R2 <- model$rows$R2
  C1 <- model$cols$C1
  C2 <- model$cols$C2
  
  need.exp <- score %in% c("t-score", "z-score")
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
    }
    else {
      scores <- (O - E) / sqrt(O + 1) # "discounted" t-score
    }
  } else if (score == "z-score") {
    if (sparse) {
      idx <- nonzero & O > E
      scores[idx] <- (O[idx] - E[idx]) / sqrt(E[idx])
    }
    else {
      scores <- (O - E) / sqrt(E)
    }
  } else if (score == "Dice") {
    scores <- 2 * O / outer(R1, C1, "+")
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

