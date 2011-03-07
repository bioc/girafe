\name{perWindow}
\Rdversion{1.1}
\alias{perWindow}
\title{Investigate aligned reads in genome intervals with sliding
  windows}
\description{
  Investigate aligned reads in genome intervals with sliding
  windows.
}
\usage{
perWindow(object, chr, winsize, step, normaliseByMatches = TRUE,
          mem.friendly = FALSE)
}
\arguments{
  \item{object}{object of class \code{AlignedGenomeIntervals}}
  \item{chr}{string; which chromosome to investigate with sliding windows}
  \item{winsize}{integer; size of the sliding window in base-pairs}
  \item{step}{integer; offset between the start positions of two sliding
    windows}
  \item{normaliseByMatches}{logical; should the number of reads per
     \code{AlignedGenomeInterval} be normalised by the number of genomic
     matches of the read sequence before summing them up in each
     window? (\emph{i.e.} derivation a weighted sum of read counts)}
  \item{mem.friendly}{logical; argument passed on to function
    \code{interval_overlap}; if \code{TRUE} the less RAM and, if the
    \code{multicore} package is attached, multiple processors are used
    for computing the overlap, on the expense of time}
}
\details{
  The windows are constructed from the first base position onto which a
  read has been mapped until the end of the chromosome.
}
\value{
  a \code{data.frame} with the following information for each sliding
  window on the chromosome
  \item{chr}{string; which chromosome the interval is on}
  \item{start}{integer; start coordinate of the windows on the
    chromosome}
  \item{end}{integer; end coordinate of the windows on the
    chromosome}
  \item{n.overlap}{integer; number of read match positions inside the
    window. Per match position there can be one or more reads mapped, so
    this number always is smaller than \code{n.reads}}
  \item{n.reads}{numeric; number of reads which match positions inside
    this window; can be floating-point numbers if argument
    \code{normaliseByMatches=TRUE}}
  \item{n.unique}{integer; number of reads which each only have one
    match position in the genome and for which this position is
    contained inside this window}
  \item{max.reads}{integer; the maximal number of reads at any single
    one match position contained inside this window}
  \item{first}{integer; coordinate of the first read alignment found
    inside the window}
  \item{last}{integer; coordinate of the last read alignment found
    inside the window}
  %\item{score}{numeric; a custom score that should indicate how relevant
  % this window is for our analyses. This score is very preliminary and
  %  not relevant for most users.}

  The result is of class \code{data.frame} and in addition of the
  (S3) class \code{slidingWindowSummary}, which may be utilized by
  follow-up functions.
}
\author{Joern Toedling}
\seealso{
  \code{\link{AlignedGenomeIntervals-class}}
}
\examples{
  exDir <- system.file("extdata", package="girafe")
  exA   <- readAligned(dirPath=exDir, type="Bowtie", 
    pattern="aravinSRNA_23_no_adapter_excerpt_mm9_unmasked.bwtmap")
  exAI  <- as(exA, "AlignedGenomeIntervals")
  exPX  <- perWindow(exAI, chr="chrX", winsize=1e5, step=0.5e5)
  head(exPX[order(exPX$n.overlap, decreasing=TRUE),])
}
\keyword{manip}
