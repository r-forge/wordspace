##
##  Collection of utility functions for illustrations
##

library(wordspace)

## Perl's qw()
qw <- function(x, sep="\\s+") unlist(strsplit(x, split=sep, perl=TRUE))

## draw arrow with start/end shortened by configurable amount and optional rotated label
draw.arrow <- function (x1, y1, x2, y2, label=NULL, label.pos=0.5, label.offset=5, cut1=0, cut2=0, head1=FALSE, head2=TRUE, length=.15, lwd=1, ...) {
  .dx <- x2 - x1
  .dy <- y2 - y1
  .l <- sqrt(.dx^2 + .dy^2)
  .code <- head2 * 1 + head1 * 2
  arrows(x1 + cut1 * .dx/.l, y1 + cut1 * .dy/.l, x2 - cut2 * .dx/.l, y2 - cut2 * .dy/.l, code=.code, length=length, lwd=lwd, ...)
  if (!missing(label)) {
    .tx <- x1 + label.pos * .dx - label.offset * .dy/.l
    .ty <- y1 + label.pos * .dy + label.offset * .dx/.l
    .angle <- 180 * atan2(.dy, .dx) / pi
    .angle <- ((.angle + 90) %% 180) - 90
    text(.tx, .ty, labels=label, srt=.angle, ...)
  }
}
