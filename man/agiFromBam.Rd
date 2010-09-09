\name{agiFromBam}
\alias{agiFromBam}
\title{
  Create AlignedGenomeIntervals objects from BAM files.
}
\description{
  Function to create \code{AlignedGenomeIntervals} objects from BAM
  (binary alignment map format) files. Uses functions from package
  \code{Rsamtools} to parse BAM files.
}
\usage{
agiFromBam(bamfile, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamfile}{
    File path of BAM file. BAM file should be sorted and have an index
    in the same directory (see Details below). }
  \item{\dots}{
    further arguments passed on to function \code{scanBam} }
}
\details{
  Note: the BAM files must be sorted and must also have an index
  file (\code{*.bai}) in the same directory. These should be done when
  creating the BAM. However, the functions \code{sortBam} and
  \code{indexBam} can be used for the same purpose, as can the
  respective modules of the \dQuote{samtools} library (\sQuote{samtools
    sort} and \sQuote{samtools index}).
  
  The BAM files are parsed chromosome by chromosome to limit the memory
  footprint of the function. Thus, this function aims to be a
  less-memory-consuming alternative to first reading in the BAM file
  using the \code{readAligned} function and then converting the
  \code{AlignedRead} object into an \code{AlignedGenomeIntervals}
  object.
}
\value{
  An object of class \code{AlignedGenomeIntervals}.
}
\references{
  \url{http://samtools.sourceforge.net}
}
\author{J Toedling}
\seealso{\code{\link[Rsamtools]{scanBam}},
  \code{\link{AlignedGenomeIntervals-class}}
}
\examples{
  fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
  ExGi <- agiFromBam(fl)
  head(detail(ExGi))
}
\keyword{manip}
