dsm.score <- function (model, score="frequency",
                       sparse=TRUE, negative.ok=NA,
                       transform=c("none", "log", "root", "sigmoid"),
                       scale=c("none", "standardize", "center", "scale"),
                       normalize=FALSE, method="euclidean", p=2,
                       matrix.only=FALSE, update.nnzero=FALSE,
                       batchsize=1e6) {
  ## validate DSM object
  model.info <- check.dsm(model, validate=TRUE)
  have.M <- model.info$M$ok
  have.S <- model.info$S$ok
  
  ## check association score and determine internal code (must be consistent with <score.c>)
  force(score) # need to evaluate argument to check whether it is a function or a string
  if (is.function(score)) {
    AM <- score     # callback function implementing user-defined AM
    score <- "user" # indicates a user-define association score
    if (!have.M) stop("cannot compute association scores: no co-occurrence frequency data available")
    if (model.info$locked) stop("marginal frequencies are invalid, cannot compute association scores")
  } else {
    score <- match.arg(score, c("frequency", "simple-ll", "t-score", "z-score", "Dice", "MI", "tf.idf", "reweight"))
    score.code <- switch(score, frequency=0, reweight=0, "simple-ll"=1, "t-score"=2, "z-score"=3, Dice=4, MI=5, tf.idf=6)
    if (score.code %in% c(0, 4, 6)) sparse <- TRUE # frequency measure, reweighting, tf.idf and Dice are always sparse
    if (score == "reweight" && !have.S) stop("cannot use score='reweight': association scores have not been computed yet")
    if (score != "reweight" && !have.M) stop("cannot compute association scores: no co-occurrence frequency data available")
    if (model.info$locked && !(score.code %in% c(0, 6))) stop("marginal frequencies are invalid, cannot compute association scores")
  }

  ## check transformation and determine internal code (must be consistent with <score.c>)
  transform <- match.arg(transform)
  transform.code <- switch(transform, none=0, log=1, root=2, sigmoid=3)
  transform.fun <- switch( # for use with user-defined AM
    transform,
    none = identity,
    log = function (x) sign(x) * log(abs(x) + 1),
    root = function (x) sign(x) * sqrt(abs(x)),
    sigmoid = tanh
  )
    
  ## check other arguments
  scale <- match.arg(scale)
 
  ## set up marginal frequencies for different AMs
  if (score == "reweight") {
    cooc.matrix <- model$S # neat trick: apply "frequency" association measure to S instead of M
    f1 <- 0 # we may need dummy entries for marginal frequencies and sample size
    f2 <- 0
    N  <- 0
  } else if (score == "tf.idf") {
    cooc.matrix <- model$M
    f1 <- model$rows$f # dummy, will be ignored
    if ("df" %in% colnames(model$cols)) {
      f2 <- model$cols$df
      N <- 1             # df should contain relative document frequencies -> dummy document count
    } else if ("nnzero" %in% colnames(model$cols)) {
      f2 <- model$cols$nnzero + 1
      N <- nrow(cooc.matrix) + 1  # simulate df by column nonzero counts, relative to number of rows of the matrix
    } else {
      stop("relative document frequencies ('df') or nonzero counts ('nnzero') for feature terms (= columns) are required in order to compute tf.idf association measure")
    } 
  } else {
    if (!have.M) stop("cannot compute association scores: no co-occurrence frequency data available")
    cooc.matrix <- model$M
    f1 <- model$rows$f
    f2 <- model$cols$f
    N  <- model.info$N
    if (is.null(f1)) stop("cannot compute association scores: no marginal frequencies for target terms")
    if (is.null(f2)) stop("cannot compute association scores: no marginal frequencies for feature terms")
    if (is.na(N)) stop("cannot compute association scores: unknown sample size")
  }

  ## check matrix format and sparse/dense representation
  matrix.info <- dsm.is.canonical(cooc.matrix)
  cooc.sparse <- matrix.info$sparse
  if (is.na(negative.ok)) negative.ok <- !cooc.sparse

  if (!negative.ok) {
    if (!sparse) stop("computation of non-sparse association scores would introduce negative values and force dense representation: specify negative.ok=TRUE if you really want to do this")
    if (scale %in% c("standardize", "center")) stop("column scaling would introduce negative values and force dense representation: specify negative.ok=TRUE if you really want to do this")
  }
  if (!sparse && cooc.sparse) {
    cooc.matrix <- as.matrix(cooc.matrix) # make co-occurrence matrix dense for non-sparse association scores
    matrix.info <- list(sparse=FALSE, canonical=TRUE, nonneg=FALSE)
    cooc.sparse <- FALSE
  }
  
  if (!matrix.info$canonical) cooc.matrix <- dsm.canonical.matrix(cooc.matrix)

  ## compute association scores and apply optional transformation
  if (score == "user") {
    ## (A) user-defined association measure: process large matrix in batches

    ## wrapper around callback function provides observed and expected frequencies as arguments with lazy evaluation
    compute.AM <- function (
      AM, f, f1, f2, N, rows, cols,
      O=f, E=f1*f2/N,
      R1=f1, R2=N-f1, C1=f2, C2=N-f2,
      O11=f, O12=f1-f, O21=f2-f, O22=N-f1-f2+f,
      E11=f1*f2/N, E12=f1*C2/N, E21=R2*f2/N, E22=R2*C2/N) {
      AM(f=f, f1=f1, f2=f2, N=N, rows=rows, cols=cols, O=O, E=E, 
         R1=R1, R2=R2, C1=C1, C2=C2,
         O11=O11, O12=O12, O21=O21, O22=O22, E11=E, E12=E12, E21=E21, E22=E22)
    }

    if (cooc.sparse) {
      ## sparse matrix: unpack dgCMatrix into triplet representation (i.row, i.col, f), then process in batches
      # i.row <- cooc.matrix@i + 1 # we compute i.row[idx] directly in the loop below
      i.col <- rep(seq_len(ncol(cooc.matrix)), times=diff(cooc.matrix@p)) # large vector, but can't be done effectively in batches
      n <- length(cooc.matrix@x)
      scores.x <- numeric(n) # pre-allocate result vector
      i1 <- 1
      while (i1 <= n) {
        i2 <- min(i1 + batchsize - 1, n)
        ## cat(sprintf("dsm.score: processing cells #%d .. #%d\n", i1, i2))
        idx <- i1:i2 # cells to process in this batch
        i.row.idx <- cooc.matrix@i[idx] + 1L
        i.col.idx <- i.col[idx]
        scores.x[idx] <- pmax(0, transform.fun(compute.AM(
          AM, f=cooc.matrix@x[idx], f1=f1[i.row.idx], f2=f2[i.col.idx], N=N,
          rows=model$rows[i.row.idx, ], cols=model$cols[i.col.idx, ]
        )))
        i1 <- i1 + batchsize
        gc() # clean up temporary objects so they don't accumulate in RAM
      }
      rm(i.col, idx); gc()
      scores <- new("dgCMatrix", Dim=as.integer(c(model.info$nrow, model.info$ncol)), p=cooc.matrix@p, i=cooc.matrix@i, x=scores.x)
      rm(scores.x); gc()
    } else {
      ## dense matrix: divide columns into batches
      nR <- nrow(cooc.matrix)
      nC <- ncol(cooc.matrix)
      batch.cols <- ceiling(batchsize / nC) # number of columns per batch
      scores <- matrix(0, nR, nC)
      i1 <- 1
      while (i1 <= nC) {
        i2 <- min(i1 + batch.cols - 1, nC)
        ## cat(sprintf("dsm.score: processing columns #%d .. #%d\n", i1, i2))
        idx <- i1:i2 # columns in batch
        i.row <- rep(1:nR, length(idx)) # row index for batch matrix
        i.col <- rep(idx, each=nR)      # column index for batch matrix
        scores[, idx] <- pmax(if (sparse) 0 else -Inf, transform.fun(compute.AM(
          AM, f=cooc.matrix[, idx, drop=FALSE], f1=f1[i.row], f2=f2[i.col], N=N, 
          rows=model$rows[i.row, ], cols=model$cols[i.col, ]
        )))
        i1 <- i1 + batch.cols
        gc() # clean up temporary objects so they don't accumulate in RAM
      }
      rm(i.row, i.col, idx); gc()
    }
  } else {
    ## (B) built-in association measures: C code for optimal memory-efficiency

    if (cooc.sparse) {
      ## compute association scores for sparse matrix (in canonical dgCMatrix format)
      scores.x <- CPP_dsm_score_sparse(model.info$nrow, model.info$ncol, cooc.matrix@p, cooc.matrix@i, cooc.matrix@x, f1, f2, N, score.code, sparse, transform.code)
      scores <- new("dgCMatrix", Dim=as.integer(c(model.info$nrow, model.info$ncol)), p=cooc.matrix@p, i=cooc.matrix@i, x=scores.x)
      rm(scores.x)
    } else {
      ## compute dense or sparse association scores for dense matrix
      scores <- CPP_dsm_score_dense(cooc.matrix, f1, f2, N, score.code, sparse, transform.code)
    }
  }

  if (scale == "standardize") {
    scores <- scale(scores, center=TRUE, scale=TRUE)
  } else if (scale == "center") {
    scores <- scale(scores, center=TRUE, scale=FALSE)
  } else if (scale == "scale") {
    rms <- colNorms(scores, "euclidean") / sqrt(nrow(scores) - 1) # root mean square according to ?scale
    scores <- scaleMargins(scores, cols = 1 / rms, duplicate=FALSE) # transform in-place (scores has been allocated above)
  } else {
    # no scaling
  }

  if (normalize) {
    ## carry out normalize.rows() operation in-place (because scores has been newly allocated above)
    scores <- scaleMargins(scores, 1 / rowNorms(scores, method=method, p=p), duplicate=FALSE)
  }

  dimnames(scores) <- dimnames(cooc.matrix) # make sure that row and column names are preserved
  if (!negative.ok) attr(scores, "nonneg") <- TRUE # S is known to be non-negative
  if (matrix.only) {
    return(scores)
  } else {
    model$S <- scores
    if (update.nnzero) {
      model$rows$nnzero <- rowNorms(scores, method="minkowski", p=0)
      model$cols$nnzero <- colNorms(scores, method="minkowski", p=0)
    }
    return(model)
  }
}

