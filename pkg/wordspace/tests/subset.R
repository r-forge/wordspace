## Validate subset() methods for DSM objects

library(wordspace)

## example from ?subset.dsm)
model <- DSM_TermContext

m2 <- subset(model, nchar(term) <= 4, nnzero <= 3)
stopifnot(all.equal(dim(m2), c(3, 4)))              # should extract a 3 x 4 matrix
stopifnot(all.equal(m2$rows$term, c("cat", "dog", "time"))) # check that we have the expected rows
stopifnot(all.equal(m2$rows$term, rownames(m2$M)))  # and that dimnames of the matrix are set correctly 
stopifnot(all.equal(m2$cols$nnzero, c(2, 2, 1, 0))) # check that nonzero counts have been updated

m2 <- subset(model, nnzero <= 3, nnzero <= 3, drop.zeroes=TRUE)
stopifnot(all.equal(dim(m2), c(3, 2))) # make sure that empty rows and columns are dropped
stopifnot(all(m2$rows$nnzero > 0) && all(m2$cols$nnzero > 0))

## without recursive=TRUE, condition on nonzero count is only satisfied in original matrix
m2 <- subset(model, nchar(term) <= 4, nnzero >= 2)
stopifnot(all.equal(dim(m2), c(3, 7)))

m2 <- subset(model, nchar(term) <= 4, nnzero >= 2, recursive=TRUE)
stopifnot(all.equal(dim(m2), c(3, 4)))
