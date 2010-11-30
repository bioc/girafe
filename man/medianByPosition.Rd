\name{medianByPosition}
\alias{medianByPosition}
\title{Compute median quality for each nucleotide position}
\description{
  This function computes the median quality for each position
  in a read over all reads in a ShortReadQ object.
}
\usage{
medianByPosition(x, method = "Sanger", batchSize = 100000L)
}
\arguments{
  \item{x}{object of class \code{ShortReadQ}, such as the result of
    function \code{readFastq}}
  \item{method}{string; passed on to function \code{intPhred}}
  \item{batchSize}{number of rows to process in each iteration; directly
    influences RAM usage of this function}
}
\details{
  The quality values are computed for each batch of reads and stored as
  numeric \code{Rle} objects for each postion. In each iteration, the
  \code{Rle} object of the current batch is merged with the previous one
  in order to keep the RAM usage low.
}
\value{
  A numeric vector of the median values per nucleotide position in the
  reads. The length of this vector corresponds to the length of the
  longest read in the data.
}
%\references{}
\author{Joern Toedling}
%\note{}
\seealso{\code{\link{intPhred}}}
\examples{
   exDir <- system.file("extdata", package="girafe")
   ra  <- readFastq(dirPath=exDir, pattern=
                   "aravinSRNA_23_plus_adapter_excerpt.fastq")
   medianByPosition(ra, batchSize=200)
}
\keyword{manip}
