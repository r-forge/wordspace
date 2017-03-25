## Test evaluation functions on real and toy tasks

library(wordspace)

qw <- function (x) unlist(strsplit(x, "\\s+", perl=TRUE)) # Perl's qw()

## a small invented multiple-choice task that works with the minimalistic DSM below
Mochi <- data.frame(
  
)

M <- dsm.score(DSM_TermTerm, score="frequency", transform="log", matrix.only=TRUE)
DM <- dist.matrix(M) # angluar distances
