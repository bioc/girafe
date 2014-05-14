\name{intPhred}
\Rdversion{1.1}
\alias{intPhred}
\title{Extract integer Phred score values from FastQ data}
\description{
  Function to extract integer Phred score values from FastQ data.
}
\usage{
intPhred(x, method="Sanger", returnType="list")
}
\arguments{
  \item{x}{object of class \code{ShortReadQ}; which contains read
    sequences and quality scores; usually read in from a Fastq files.}
  \item{method}{string; one of 'Sanger', 'Solexa' or
    'previousSolexa'. See details below.}
  \item{returnType}{string; in which format should the result
    be returned, either as a 'list' or as a 'matrix'.}
}
\details{
  There are different standards for encoding read qualities in Fastq
  files. The 'Sanger' format encodes a Phred quality score from 0 to 93
  using ASCII 33 to 126. The current 'Solexa'/llumina format (1.3 and
  higher) encodes a Phred quality score from 0 to 40 using ASCII 64 to
  104. The 'previous Solexa'/Illumina format (1.0) encodes a custom
  Solexa/Illumina quality score from -5 to 40 using ASCII 59 to 104.
  This custom Solexa quality score is approximately equal to the Phred
  scores for high qualities, but differs in the low quality range.
}
\value{
  If \code{returnType} is equal to \sQuote{list}:
  A list of integer Phred quality values of the same length as the
  number of reads in the object \code{x}.
  
  If \code{returnType} is equal to \sQuote{matrix}:
  A matrix of integer Phred quality values. The number of rows is the
  number of reads in the object \code{x}. The number of columns is the
  maximum length (width) over all reads in object \code{x}. The last
  entries for reads that are shorter than this maximum width are
  'NA'.
}
\references{
  \url{http://maq.sourceforge.net/fastq.shtml}
}
\author{Joern Toedling}
\seealso{\code{\link[ShortRead]{ShortReadQ-class}},
  \code{\link[ShortRead]{readFastq}}}
\examples{
  exDir <- system.file("extdata", package="girafe")
  ra  <- readFastq(dirPath=exDir, pattern=
            "aravinSRNA_23_plus_adapter_excerpt.fastq")
  ra.quals <- intPhred(ra, method="Sanger",
                       returnType="matrix")
  ra.qmed  <- apply(ra.quals, 2, median)
  if (interactive())
     plot(ra.qmed, type="h", ylim=c(0,42), xlab="Base postion",
          ylab="Median Phred Quality Score", lwd=2, col="steelblue")
}
\keyword{manip}
