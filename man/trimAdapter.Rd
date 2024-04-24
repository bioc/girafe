\name{trimAdapter}
\Rdversion{1.1}
\alias{trimAdapter}
\title{Remove 3' adapter contamination}
\description{
  Function to remove 3' adapter contamination from reads
}
\usage{
trimAdapter(fq, adapter, match.score = 1, mismatch.score = -1,
            score.threshold = 2)
}
\arguments{
  \item{fq}{Object of class \code{ShortReadQ}; the reads with possible
    adapter contamination.}
  \item{adapter}{object of class \code{DNAString} or class
    \code{character}; the sequence of the 3' adapter which could give
    rise to the 3' contamination. If of class \code{character}, it is
    converted to a \code{DNAString} inside the function.}
  \item{match.score}{numeric; alignment score for matching bases}
  \item{mismatch.score}{numeric; alignment score for mismatches}
  \item{score.threshold}{numeric; minimum total alignment score required
    for an overlap match between the 3' end of the read and the 5' end
    of the adapter sequence.}
}
\details{
  Performs an overlap alignment between the ends of the reads and the
  start of the adapter sequence.
}
\note{
  The function \code{trimLRPatterns} from package \code{ShortRead} may
  be a faster alternative to this function.
}
\value{
  An object of class \code{ShortReadQ} containing the reads without the
  3' adapter contamination.
}
\author{J. Toedling}
\seealso{\code{\link[pwalign]{pairwiseAlignment}},
  \code{\link[IRanges]{narrow}},
  \code{\link[ShortRead]{readFastq}},
  \code{\link[ShortRead]{writeFastq}}
}
\examples{
  exDir <- system.file("extdata", package="girafe")
  ## load reads containing adapter fragments at the end
  ra23.wa  <- readFastq(dirPath=exDir, pattern=
              "aravinSRNA_23_plus_adapter_excerpt.fastq")
  table(width(ra23.wa))
  # adapter sequence obtained from GEO page
  #  accession number: GSE10364
  #adapter <- DNAString("CTGTAGGCACCATCAAT")
  adapter <- "CTGTAGGCACCATCAAT"

  # trim adapter
  ra23.na  <- trimAdapter(ra23.wa, adapter)
  table(width(ra23.na))
}
\keyword{manip}
