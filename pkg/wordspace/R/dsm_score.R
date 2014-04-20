dsm.score <- function (model,
                       score=c("frequency", "simple-ll", "t-score", "z-score", "Dice", "MI", "tf.idf", "reweight"),
                       sparse=TRUE, negative.ok=NA,
                       transform=c("none", "log", "root", "sigmoid"),
                       scale=c("none", "standardize", "center", "scale"),
                       normalize=FALSE, method="euclidean", p=2,
                       matrix.only=FALSE, update.nnzero=FALSE) {
  score <- match.arg(score)
  transform <- match.arg(transform)
  scale <- match.arg(scale)
  model.info <- check.dsm(model, validate=TRUE)
  have.M <- model.info$M$ok
  have.S <- model.info$S$ok
  calculate.AM <- !(score %in% c("frequency", "reweight", "tf.idf"))
  if (model.info$locked && calculate.AM) stop("marginal frequencies are invalid, cannot compute association scores")
  
  ## internal codes for association scores and transformation functions (must match C code in <score.c>)
  score.code <- switch(score, frequency=0, reweight=0, "simple-ll"=1, "t-score"=2, "z-score"=3, Dice=4, MI=5, tf.idf=6)
  transform.code <- switch(transform, none=0, log=1, root=2, sigmoid=3)
  if (score.code %in% c(0, 4, 6)) sparse <- TRUE # frequency measure, reweighting, tf.idf and Dice are always sparse

  if (score == "reweight") {
    if (!have.S) stop("cannot use score='reweight': association scores have not been computed yet")
    cooc.matrix <- model$S # neat trick: apply "frequency" association measure to S instead of M
    f1 <- 0 # we may need dummy entries for marginal frequencies and sample size
    f2 <- 0
    N  <- 0
  } else if (score == "tf.idf") {
    if (!have.M) stop("cannot compute association scores: no co-occurrence frequency data available")
    cooc.matrix <- model$M
    f1 <- model$rows$f # dummy, will be ignored
    if ("df" %in% colnames(model$cols)) {
      f2 <- model$cols$df
      N <- 1             # df should contain relative document frequencies -> dummy document count
    } else if ("nnzero" %in% colnames(model$cols)) {
      f2 <- model$cols$nnzero
      N <- nrow(cooc.matrix) # simulate df by column nonzero counts, relative to number of rows of the matrix
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

  ## -- compute association scores and apply optional transformation (C code for optimal memory-efficiency) --
  if (cooc.sparse) {
    ## compute association scores for sparse matrix (in canonical dgCMatrix format)
    scores.x <- double(length(cooc.matrix@x))
    scores.x <- CPP_dsm_score_sparse(model.info$nrow, model.info$ncol, cooc.matrix@p, cooc.matrix@i, cooc.matrix@x, f1, f2, N, score.code, sparse, transform.code)
    scores <- new("dgCMatrix", Dim=as.integer(c(model.info$nrow, model.info$ncol)), p=cooc.matrix@p, i=cooc.matrix@i, x=scores.x)
    rm(scores.x)
  } else {
    ## compute dense or sparse association scores for dense matrix
    scores <- CPP_dsm_score_dense(cooc.matrix, f1, f2, N, score.code, sparse, transform.code)
  }

  if (scale == "standardize") {
    scores <- scale(scores, center=TRUE, scale=TRUE)
  } else if (scale == "center") {
    scores <- scale(scores, center=TRUE, scale=FALSE)
  } else if (scale == "scale") {
    scores <- scaleMargins(scores, cols = 1 / colNorms(scores, "euclidean")) # preserves sparse matrix representation
  } else {
    # no scaling
  }

  if (normalize) {
    scores <- normalize.rows(scores, method=method, p=p)
  }

  dimnames(scores) <- dimnames(cooc.matrix) # make sure that row and column names are preserved
  if (matrix.only) {
    return(scores)
  } else {
    model$S <- scores
    if (update.nnzero) {
      is.nzero <- scores != 0
      model$rows$nnzero <- rowSums(is.nzero)
      model$cols$nnzero <- colSums(is.nzero)
      rm(is.nzero)
    }
    return(model)
  }
}

