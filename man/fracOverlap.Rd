\name{fracOverlap}
\Rdversion{1.1}
\alias{fracOverlap}
\title{Retrieve intervals overlapping by fraction of width}
\description{
  Function to retrieve overlapping intervals that overlap at least by a
  specified fraction of their widths.
}
\usage{fracOverlap(I1, I2, min.frac = 0.66)}
\arguments{
  \item{I1}{object that inherits from class \code{Genome_intervals}}
  \item{I2}{object that inherits from class \code{Genome_intervals}}
  \item{min.frac}{numeric; minimum fraction of smaller of the two
    interval widths by which two intervals should overlap in order to be
    marked as overlapping.}
}
\value{
  An object of class \code{data.frame} with one row each for a pair of
  overlapping elements.
  \item{Index1}{Index of interval in first interval list}
  \item{Index2}{Index of interval in second interval list}
  \item{n}{number of bases that the two intervals overlap}
  \item{fraction}{fraction of the smaller of the two intervals' width by
    which the two intervals overlap}
  \item{frac1}{fraction of interval 1's width by
    which the two intervals overlap}
  \item{frac2}{fraction of interval 2's width by
    which the two intervals overlap}
}
\author{J. Toedling}
\seealso{\code{\link[genomeIntervals]{interval_overlap}}}
\examples{
   data("gen_ints", package="genomeIntervals")
   i[4,2] <- 13
   fracOverlap(i, i)
}
\keyword{manip}
