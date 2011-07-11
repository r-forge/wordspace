dsm.score <- function (model,
                       score=c("frequency", "simple-ll", "t-score", "z-score", "Dice", "MI", "reweight"),
                       sparse=FALSE,
                       transform=c("none", "log", "root", "sigmoid"),
                       scale=c("none", "standardize", "center"),
                       normalize=FALSE, method="euclidean", p=2,
                       matrix.only=FALSE) {
  score <- match.arg(score)
  transform <- match.arg(transform)
  scale <- match.arg(scale)
  model.info <- check.dsm(model, check.S=(score=="reweight"))
  sparse.M <- model.info$sparse
  calculate.AM <- !(score %in% c("frequency", "reweight"))

  ## internal codes for scoring measures and transformation functions (must match C code in <score.c>)
  score.code <- switch(score, frequency=0, "simple-ll"=1, "t-score"=2, "z-score"=3, Dice=4, MI=5)
  transform.code <- switch(transform, none=0, log=1, root=2, sigmoid=3)
  if (score == "frequency") sparse <- TRUE # frequency is always a sparse score

  if (model.info$locked && calculate.AM) stop("marginal frequencies have been adjusted, cannot recompute association scores")
  if (score == "reweight" && !model.info$have.S) stop("cannot use score='reweight': association scores have not been computed yet")
  if (scale != "none" && sparse) warning("scaling of matrix columns destroys sparse representation: are you sure you want to do this?")
  if (normalize && !sparse) warning("normalization of matrix rows is only sensible for a sparse non-negative representation")
  if (sparse.M && calculate.AM && !sparse) stop("association scores for sparse DSM can only be computed with sparse=TRUE")
  
  if (score == "reweight") {
    scores <- model$S # reweight existing scores (with transformation, scaling and/or normalisation)
  } else if (sparse.M) {
    ## check.dsm() has already verified that model$M is of class dgCMatrix
    scores.x <- double(length(model$M@x))
    .C(
      C_dsm_score_sparse,
      scores.x,
      as.integer(model.info$ncol),
      as.integer(model$M@p),
      as.integer(model$M@i),
      as.double(model$M@x),
      as.double(model$rows$f),
      as.double(model$cols$f),
      as.double(model.info$N),
      as.integer(score.code),
      as.logical(sparse),
      as.integer(transform.code),
      DUP=FALSE, NAOK=FALSE
    )
    scores <- new("dgCMatrix", Dim=as.integer(c(model.info$nrow, model.info$ncol)), p=model$M@p, i=model$M@i, x=scores.x)
    rm(scores.x)
  } else if (!sparse.M) {
    scores <- double(model.info$nrow * model.info$ncol)
    .C(
      C_dsm_score_dense,
      scores,
      as.integer(model.info$nrow),
      as.integer(model.info$ncol),
      as.double(model$M),
      as.double(model$rows$f),
      as.double(model$cols$f),
      as.double(model.info$N),
      as.integer(score.code),
      as.logical(sparse),
      as.integer(transform.code),
      DUP=FALSE, NAOK=FALSE
    )
    dim(scores) <- c(model.info$nrow, model.info$ncol)
    ## -- remaining clauses are effectively disabled now --
  } else if (score == "frequency") {
  scores <- model$M # "frequency" AM = observed co-occurrence frequency
  }
  else {
    # need to compute association scores, based on marginals and/or expected frequencies
    
    if (sparse) {
      # -- operate on sparse representation (even if original matrix is dense)
      if (sparse.M) {
        M <- as(model$M, "dgTMatrix")
      } else {
        idx <- which(model$M > 0, arr.ind=TRUE)
        M <- new("dgTMatrix", Dim=dim(model$M), i=idx[,1]-1L, j=idx[,2]-1L, x=model$M[idx])
        rm(idx)
      }

      O <-M@x                       # observed frequencies (notation follows Evert 2004, 2008)
      N <- model$N                  # sample size
      R1 <- model$rows$f[M@i + 1]   # row marginals (compute R2 = N - R1 only if needed)
      C1 <- model$cols$f[M@j + 1]   # column marginals (compute C2 = N - C1 only if needed)
      # expected frequencies will be calculated when necessary

      if (score == "simple-ll") {
        E <- R1 * C1 / N
        score.vec <- 2 * (O * log(O / E) - (O - E))
        score.vec[O < E] <- 0
        rm(E)
      } else if (score == "t-score") {
        E <- R1 * C1 / N
        score.vec <- (O - E) / sqrt(O)
        score.vec[O < E] <- 0
        rm(E)
      } else if (score == "z-score") {
        E <- R1 * C1 / N
        score.vec <- (O - E) / sqrt(E)
        score.vec[O < E] <- 0
        rm(E)
      } else if (score == "Dice") {
        score.vec <- 2 * O / (R1 + C1)
      } else if (score == "MI") {
        E <- R1 * C1 / N
        score.vec <- log2(O / E)
        score.vec[O < E] <- 0
        rm(E)
      }
      else stop("internal error - sparse AM not implemented")
      rm(O, R1, C1)

      M@x <- score.vec
      scores <- if (sparse.M) as(M, "dgCMatrix") else as.matrix(M)
      rm(M, score.vec)
    }
    else {
      # -- standard dense matrix processing
      stopifnot(!sparse.M)

      O <- model$M                  # observed frequencies (notation follows Evert 2004, 2008)
      N <- model.info$N             # sample size
      R1 <- model$rows$f            # row marginals (compute R2 = N - R1 only if needed)
      C1 <- model$cols$f            # column marginals (compute C2 = N - C1 only if needed)

      need.exp <- score %in% c("simple-ll", "t-score", "z-score", "MI")
      if (need.exp) E <- outer(R1, C1) / N

      if (score == "simple-ll") {
        OmE <- O - E
        scores <- sign(OmE) * 2 * (ifelse(O > 0, O * log(O/E), 0) - OmE)
      } else if (score == "t-score") {
        scores <- (O - E) / sqrt(O + 1) # "discounted" t-score
      }
      else if (score == "z-score") {
        scores <- (O - E) / sqrt(E)  # z-score (E=0 should never happen)
      }
      else if (score == "Dice") {
        scores <- 2 * O / outer(R1, C1, "+") # Dice coefficient
      } else if (score == "MI") {
        scores <- log2(O / E)        # standard pointwise MI (not clear how to avoid -Inf here)
      }
      
      rm(O, R1, C1)
      if (need.exp) rm(E)
    }

  }

  # transformation and scaling now operates in the same way on dense or sparse matrix "scores"
  if (score == "reweight") {
    ## -- otherwise, transformation has already been applied by C code --
    if (transform == "log") {
      if (!sparse.M) {
        scores <- sign(scores) * log(abs(scores) + 1) # "signed log" transformation
      } else {
        scores@x <- sign(scores@x) * log(abs(scores@x) + 1)
      }
    } else if (transform == "root") {
      scores <- sign(scores) * sqrt(abs(scores)) # signed square root
    } else if (transform == "sigmoid") {
      scores <- tanh(scores) # tanh sigmoid function, approximates signed binary vector {-1, 0, +1}
    } else {
      # no transformation
    }
  }

  if (scale == "standardize") {
    scores <- scale(scores, center=TRUE, scale=TRUE)
  } else if (scale == "center") {
    scores <- scale(scores, center=TRUE, scale=FALSE)
  } else {
    # no scaling
  }

  if (normalize) {
    scores <- normalize.rows(scores, method=method, p=p)
  }

  dimnames(scores) <- dimnames(model$M) # make sure that row and column names are preserved
  if (matrix.only) {
    return(scores)
  } else {
    model$S <- scores
    return(model)
  }
}

## Old pure-R version for comparison / benchmarking purposes
dsm.score.old <- function (model,
                       score=c("frequency", "simple-ll", "t-score", "z-score", "Dice", "MI", "reweight"),
                       sparse=FALSE,
                       transform=c("none", "log", "root", "sigmoid"),
                       scale=c("none", "standardize", "center"),
                       normalize=FALSE, method="euclidean", p=2,
                       matrix.only=FALSE) {
  score <- match.arg(score)
  transform <- match.arg(transform)
  scale <- match.arg(scale)
  model.info <- check.dsm(model)
  sparse.M <- inherits(model$M, "Matrix")
  calculate.AM <- !(score %in% c("frequency", "reweight"))

  if (model.info$locked && calculate.AM) stop("marginal frequencies have been adjusted, cannot recompute association scores")
  if (score == "reweight" && !model.info$have.S) stop("cannot use score='reweight': association scores have not been computed yet")
  if (scale != "none" && sparse) warning("scaling of matrix columns destroys sparse representation: are you sure you want to do this?")
  if (normalize && !sparse) warning("normalization of matrix rows is only sensible for a sparse non-negative representation")
  if (sparse.M && calculate.AM && !sparse) stop("association scores for sparse DSM can only be computed with sparse=TRUE")
  
  if (score == "reweight") {
    scores <- model$S # reweight existing scores (with transformation, scaling and/or normalisation)
  } else if (score == "frequency") {
    scores <- model$M # "frequency" AM = observed co-occurrence frequency
  }
  else {
    # need to compute association scores, based on marginals and/or expected frequencies
    
    if (sparse) {
      # -- operate on sparse representation (even if original matrix is dense)
      if (sparse.M) {
        M <- as(model$M, "dgTMatrix")
      } else {
        idx <- which(model$M > 0, arr.ind=TRUE)
        M <- new("dgTMatrix", Dim=dim(model$M), i=idx[,1]-1L, j=idx[,2]-1L, x=model$M[idx])
        rm(idx)
      }

      O <-M@x                       # observed frequencies (notation follows Evert 2004, 2008)
      N <- model$N                  # sample size
      R1 <- model$rows$f[M@i + 1]   # row marginals (compute R2 = N - R1 only if needed)
      C1 <- model$cols$f[M@j + 1]   # column marginals (compute C2 = N - C1 only if needed)
      # expected frequencies will be calculated when necessary

      if (score == "simple-ll") {
        E <- R1 * C1 / N
        score.vec <- 2 * (O * log(O / E) - (O - E))
        score.vec[O < E] <- 0
        rm(E)
      } else if (score == "t-score") {
        E <- R1 * C1 / N
        score.vec <- (O - E) / sqrt(O)
        score.vec[O < E] <- 0
        rm(E)
      } else if (score == "z-score") {
        E <- R1 * C1 / N
        score.vec <- (O - E) / sqrt(E)
        score.vec[O < E] <- 0
        rm(E)
      } else if (score == "Dice") {
        score.vec <- 2 * O / (R1 + C1)
      } else if (score == "MI") {
        E <- R1 * C1 / N
        score.vec <- log2(O / E)
        score.vec[O < E] <- 0
        rm(E)
      }
      else stop("internal error - sparse AM not implemented")
      rm(O, R1, C1)

      M@x <- score.vec
      scores <- if (sparse.M) as(M, "dgCMatrix") else as.matrix(M)
      rm(M, score.vec)
    }
    else {
      # -- standard dense matrix processing
      stopifnot(!sparse.M)

      O <- model$M                  # observed frequencies (notation follows Evert 2004, 2008)
      N <- model.info$N             # sample size
      R1 <- model$rows$f            # row marginals (compute R2 = N - R1 only if needed)
      C1 <- model$cols$f            # column marginals (compute C2 = N - C1 only if needed)

      need.exp <- score %in% c("simple-ll", "t-score", "z-score", "MI")
      if (need.exp) E <- outer(R1, C1) / N

      if (score == "simple-ll") {
        OmE <- O - E
        scores <- sign(OmE) * 2 * (ifelse(O > 0, O * log(O/E), 0) - OmE)
      } else if (score == "t-score") {
        scores <- (O - E) / sqrt(O + 1) # "discounted" t-score
      }
      else if (score == "z-score") {
        scores <- (O - E) / sqrt(E)  # z-score (E=0 should never happen)
      }
      else if (score == "Dice") {
        scores <- 2 * O / outer(R1, C1, "+") # Dice coefficient
      } else if (score == "MI") {
        scores <- log2(O / E)        # standard pointwise MI (not clear how to avoid -Inf here)
      }
      
      rm(O, R1, C1)
      if (need.exp) rm(E)
    }

  }

  # transformation and scaling now operates in the same way on dense or sparse matrix "scores"
  if (transform == "log") {
    if (!sparse.M) {
      scores <- sign(scores) * log(abs(scores) + 1) # "signed log" transformation
    } else {
      scores@x <- sign(scores@x) * log(abs(scores@x) + 1)
    }
  } else if (transform == "root") {
    scores <- sign(scores) * sqrt(abs(scores)) # signed square root
  } else if (transform == "sigmoid") {
    scores <- tanh(scores) # tanh sigmoid function, approximates signed binary vector {-1, 0, +1}
  } else {
    # no transformation
  }

  if (scale == "standardize") {
    scores <- scale(scores, center=TRUE, scale=TRUE)
  } else if (scale == "center") {
    scores <- scale(scores, center=TRUE, scale=FALSE)
  } else {
    # no scaling
  }

  if (normalize) {
    scores <- normalize.rows.old(scores, method=method, p=p)
  }

  dimnames(scores) <- dimnames(model$M) # make sure that row and column names are preserved
  if (matrix.only) {
    return(scores)
  } else {
    model$S <- scores
    return(model)
  }
}

