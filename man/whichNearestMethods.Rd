\name{which_nearest-methods}
\docType{methods}
\alias{which_nearest}
\alias{which_nearest,AlignedGenomeIntervals,Genome_intervals_stranded-method}
\alias{which_nearest,Genome_intervals,Genome_intervals-method}
\alias{which_nearest,Genome_intervals_stranded,Genome_intervals_stranded-method}
\title{Methods for function 'which_nearest' and genome intervals}
\description{
  For each genome interval in one set, finds the nearest interval in a
  second set of genome intervals.
}
\section{Methods}{
  \describe{
    Currently, the package \emph{girafe} contains method implementations
    for the first object (Argument: \code{from}) being of any of the
    classes \dQuote{AlignedGenomeIntervals},\dQuote{Genome_intervals} or
    \dQuote{Genome_intervals_stranded}. The second object  (Argument:
    \code{to}) has be of class \dQuote{Genome_intervals_stranded} or
    \dQuote{Genome_intervals}.
}}
\value{
  a \code{data.frame} with a number of rows equal to the number of
  intervals in argument \code{from}. The elements of the data.frame are:
  \item{distance_to_nearest}{numeric; distance to nearest interval from
    object \code{to}. Is 0 if the current interval in object \code{from}
    did overlap one or more intervals in object \code{to}}
  \item{which_nearest}{list; each list element are the indices or the
    intervals in object \code{to} that have the closest distance to the
    current interval in object \code{from}}
 \item{which_overlap}{list; each list element are the indices or the
    intervals in object \code{to} that do overlap with the
    current interval in object \code{from}}
}
\note{
  If the supplied objects are stranded, as it is the case with objects
  of classes \sQuote{AlignedGenomeIntervals} and
  \sQuote{Genome_intervals_stranded}, then the overlap and distance is
  solely computed between intervals on the same strand.

  For objects of class \sQuote{Genome_intervals}, overlap and distances
  are computed regardless of strand information.
}
\author{Joern Toedling}
\seealso{\code{\link[intervals]{which_nearest}}}
\examples{
  ### process aligned reads
  exDir <- system.file("extdata", package="girafe")
  exA   <- readAligned(dirPath=exDir, type="Bowtie", 
    pattern="aravinSRNA_23_no_adapter_excerpt_mm9_unmasked.bwtmap")
  exAI <- as(exA, "AlignedGenomeIntervals")

  ## load annotated genome features
  load(file.path(exDir, "mgi_gi.RData"))

  ## subset for sake of speed:
  A <- exAI[is.element(seqnames(exAI), c("chrX","chrY"))]
  G <- mgi.gi[is.element(seqnames(mgi.gi), c("chrX","chrY"))]

  ## find nearest annotated feature for each AlignedGenomeInterval
  WN <- which_nearest(A, G)
  dim(WN); tail(WN)

  ## notice the difference to:
  tail(which_nearest(as(A, "Genome_intervals"), G))
  # the last interval in A is located antisense to a gene,
  # but not overlapping anything on the same strand
}
\keyword{methods}
