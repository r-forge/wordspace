## Test various convenience methods for DSM objects and related functionality

library(wordspace)

dsm1 <- DSM_TermContext
dsm2 <- subset(DSM_TermTerm, f < 100000) # reduce to 5x7 matrix
dsm2 <- dsm.score(dsm2, function(O, E, ...) round(log2(O/E), 2))


## dim() of DSM object
stopifnot(all.equal(dim(dsm1), c(7, 7)))
stopifnot(all.equal(nrow(dsm2), 5))
stopifnot(all.equal(ncol(dsm2), 7))


## reading and setting dimnames
stopifnot(all.equal(dimnames(dsm1), dimnames(dsm1$M)))
stopifnot(all.equal(dimnames(dsm2), dimnames(dsm2$M)))
rn <- c("cat", "dog", "animal", "reason", "cause") # verify row/column names directly
cn <- c("breed", "tail", "feed", "kill", "important", "explain", "likely")
stopifnot(all.equal(rownames(dsm2), rn))
stopifnot(all.equal(colnames(dsm2), cn))

dsm2n <- dsm2 # modification of row/column names
rownames(dsm2n)[3] <- "pet"
stopifnot(all.equal(dimnames(dsm2), list(rn, cn))) # must not affect original object
stopifnot(all.equal(colnames(dsm2n), cn)) # no changes to column names yet
rn2 <- rn; rn2[3] <- "pet" 
stopifnot(all.equal(rownames(dsm2n), rn2)) # but one element of rownames has been changed
cn2 <- LETTERS[1:7] # replace all column names
colnames(dsm2n) <- cn2
stopifnot(all.equal(dimnames(dsm2), list(rn, cn)))
stopifnot(all.equal(rownames(dsm2n), rn2))
stopifnot(all.equal(colnames(dsm2n), cn2))

invisible(check.dsm(dsm2n, validate=TRUE)) # checks that rows$term and cols$term have been updated


## extraction of co-occurrence or score matrix
stopifnot(all.equal(as.matrix(dsm1), dsm1$M))
stopifnot(all.equal(as.matrix(dsm2), dsm2$S)) # automatic selection
stopifnot(all.equal(as.matrix(dsm2, what="M"), dsm2$M))
stopifnot(all.equal(as.matrix(dsm2, what="S"), dsm2$S))


## transposition
dsm2t <- t(dsm2)
stopifnot(nrow(dsm2t) == ncol(dsm2))
stopifnot(ncol(dsm2t) == nrow(dsm2))
stopifnot(all.equal(dsm2t$M, t(dsm2$M)))
stopifnot(all.equal(dsm2t$S, t(dsm2$S)))
