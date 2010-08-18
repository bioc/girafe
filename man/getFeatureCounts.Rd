\name{getFeatureCounts}
\alias{getFeatureCounts}
\title{get the read counts for a supplied set of genomic features}
\description{
  get the read counts for a supplied set of genomic features
}
\usage{
getFeatureCounts(AI, FG, nameColumn = "Name", fractionIncluded = 1,
returnType = "AlignedGenomeIntervals")
}
\arguments{
  \item{AI}{\code{AlignedGenomeIntervals} object}
  \item{FG}{\code{Genome\_intervals} objects of genomic features}
  \item{nameColumn}{character indicating which column of the object
    \code{FG} holds the identifiers of the genomic features; is used to
    assess the number of genomic copies per feature}
  \item{fractionIncluded}{double; which fraction of an interval needs to
    be included in a feature in order to count for the feature}
  \item{returnType}{one of \code{AlignedGenomeIntervals} or
    \code{integer}}
}
%\details{}
\value{
  depends on argument \code{returnType}:
  one of \code{AlignedGenomeIntervals} or
  a named \code{integer}
}
\author{Joern Toedling}
\seealso{\code{\link{fracOverlap}}}
%\examples{}
\keyword{manip}
\keyword{internal}% __ONLY ONE__ keyword per line
