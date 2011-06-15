# define .onLoad and .onAttach here if package initialisation functions are needed;
# .Last.lib (has to be exported) or .onUnload for package finalisation

R_sqrt <- function (x) {
  .x <- as.double(x)
  .l <- length(x)
  .res <- double(.l)
  .C(C_sqrt, .x, .l, result=.res, DUP=FALSE, NAOK=FALSE)$result
}