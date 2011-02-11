\name{weightedConsensusMatrix}
\alias{weightedConsensusMatrix}
\title{compute weighted consensus matrix}
\description{
  computes weighted consensus matrix
}
\usage{
weightedConsensusMatrix(seqs, weights, shift = NULL,
                        baseLetters = c("A", "C", "G", "T", "N"))
}
\arguments{
  \item{seqs}{character vector of read sequences}
  \item{weights}{integer; weights (read counts)}
  \item{shift}{integer; shift of each read sequence relative to the
    first column of the consensus matrix, by default: 0}
  \item{baseLetters}{alphabet}
}
\value{
  A consensus matrix
}
\author{J Toedling}

\seealso{\code{\link[Biostrings]{consensusMatrix}}}
\examples{
 ## Align following sequences with weights:
 ##   ACATT       1
 ##    CGTTA     10
 ##     TTG       3
 ##  GACATT       4

 dweights <- c(1L, 10L, 3L, 4L)
 d <- c("ACATT","CGTTA", "TTG", "GACATT")
 dshifts <- c(0L, 1L, 2L, -1L)
 W <- weightedConsensusMatrix(d, dweights, shift=dshifts)
 consensusString(W, ambiguityMap="N")
}
\keyword{manip}
\keyword{internal}
