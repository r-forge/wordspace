##
##  Artificial example matrices (defined inline)
##  .rda files are directly written to package data/ directory
##

library(wordspace)

##
## Hieroglyphs example (with English labels)
##
hieroglyphs.txt <- "get see use hear eat kill
  knife   51  20  84    0   3    0
  cat     52  58   4    4   6   26
  dog    115  83  10   42  33   17
  boat    59  39  23    4   0    0
  cup     98  14   6    2   1    0
  pig     12  17   3    2   9   27
  banana  11   2   2    0  18    0
"
fh <- textConnection(hieroglyphs.txt)
DSM_HieroglyphsMatrix <- as.matrix(read.table(fh))
mode(DSM_HieroglyphsMatrix) <- "double"
close(fh)

print(DSM_HieroglyphsMatrix)
str(DSM_HieroglyphsMatrix)
save(DSM_HieroglyphsMatrix, file="../pkg/wordspace/data/DSM_HieroglyphsMatrix.rda", compress="gzip")


##
## Wikipedia term-context matrix
##
tc.txt <- "Felidae Pet Feral Boat Philosophy Kant Back_Pain
cat             10  10     7    0          0    0         0
dog              0  10     4   11          0    0         0
animal           2  15    10    2          0    0         0
time             1   0     0    0          2    1         0
reason           0   1     0    0          1    4         1
cause            0   0     0    2          1    2         6
effect           0   0     0    1          0    1         0
"
fh <- textConnection(tc.txt)
DSM_TermContextMatrix <- as.matrix(read.table(fh))
mode(DSM_TermContextMatrix) <- "double"
close(fh)

print(DSM_TermContextMatrix)
str(DSM_TermContextMatrix)
save(DSM_TermContextMatrix, file="../pkg/wordspace/data/DSM_TermContextMatrix.rda", compress="gzip")


##
## Wikipedia term-term matrix
##
tt.txt <- "breed  tail  feed  kill  important  explain  likely
cat           83    17     7    37          0        1       0
dog          561    13    30    60          1        2       4
animal        42    10   109   134         13        5       5
time          19     9    29   117         81       34     109
reason         1     0     2    14         68      140      47
cause          0     1     0     4         55       34      55
effect         0     0     1     6         60       35      17
"
fh <- textConnection(tt.txt)
DSM_TermTermMatrix <- as.matrix(read.table(fh))
mode(DSM_TermTermMatrix) <- "double"
close(fh)

print(DSM_TermTermMatrix)
str(DSM_TermTermMatrix)
save(DSM_TermTermMatrix, file="../pkg/wordspace/data/DSM_TermTermMatrix.rda", compress="gzip")
