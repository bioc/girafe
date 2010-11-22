\name{fracOverlap}
\Rdversion{1.1}
\alias{fracOverlap}
\title{Retrieve intervals overlapping by fraction of width}
\description{
  Function to retrieve overlapping intervals that overlap at least by a
  specified fraction of their widths.
}
\usage{fracOverlap(I1, I2, min.frac=0.0, both=TRUE, mem.friendly=FALSE)}
\arguments{
  \item{I1}{object that inherits from class \code{Genome_intervals}}
  \item{I2}{object that inherits from class \code{Genome_intervals}}
  \item{min.frac}{numeric; minimum required fraction of each of the two
    interval widths by which two intervals should overlap in order to be
    marked as overlapping.}
  \item{both}{logical; shall both overlap partners meet the minimum
    fraction {\code{min.frac}} requirement? If \code{FALSE}, then
    overlaps with only partner involved to at least that fraction are
    also reported.}
  \item{mem.friendly}{logical; if set to \code{TRUE} an older but
    memory-friendlier version of \code{interval_overlap} is used
    inside this function. Note that \code{mem.friendly} is only
    evaluated if \code{I1} or \code{I2} is of class
    \code{AlignedGenomeIntervals}.}
}
\value{
  An object of class \code{data.frame} with one row each for a pair of
  overlapping elements.
  \item{Index1}{Index of interval in first interval list}
  \item{Index2}{Index of interval in second interval list}
  \item{n}{number of bases that the two intervals overlap}
  \item{fraction1}{fraction of interval 1's width by
    which the two intervals overlap}
  \item{fraction2}{fraction of interval 2's width by
    which the two intervals overlap}
}
\author{J. Toedling}
\seealso{\code{\link[genomeIntervals]{interval_overlap}}}
\examples{
   data("gen_ints", package="genomeIntervals")
   i[4,2] <- 13L
   fracOverlap(i, i, 0.5)
}
\keyword{manip}
