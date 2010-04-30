\name{plotAligned}
\Rdversion{1.1}
\alias{plotAligned}
\title{Visualize reads aligned to genome intervals}
\description{
  Visualize reads aligned to genome intervals
}
\usage{
plotAligned(x, y, chr, start, end, plus.col = "#00441b",
   minus.col = "#283d78", gff, featureLegend = FALSE,
   gffChrColumn = "seq_id", gffNameColumn="name",
   featureExclude = c("chromosome", "nucleotide_match", "insertion"),
   show="both", ylim, highlight, main, ...)
}
\arguments{
  \item{x}{Object of class \code{AlignedGenomeIntervals}}
  \item{y}{This argument is only specified for compability with
    plot.default and not used in the function.}
  \item{chr}{string; on which chromosome is the region to plot}
  \item{start}{integer; start coordinate of the chromosome region to
    plot} 
  \item{end}{integer; end coordinate of the chromosome region to
    plot} 
  \item{plus.col}{
    which colour to use for the reads on the Plus strand}
  \item{minus.col}{
     which colour to use for the reads on the Plus strand}
\item{gff}{Data frame containing annotation for genomic feature to be
    used to further annotate the plot.}
  \item{featureLegend}{logical; should a legend that describes the
    colour code for the annotated genome features be appended at the
    bottom of the plot?}
  \item{gffChrColumn}{string; which column of the \code{gff} data.frame
    holds the chromosome identifier of each feature.}
  \item{gffNameColumn}{what is the column of the gff data.frame called
    that holds the identifier of the element that should be displayed in
    the plot; default: \dQuote{name}}
  \item{featureExclude}{character; which kinds of annotated genome
    features specified in the \code{gff} are to be ignored for the plot}
  \item{show}{string; which strands to show in the plot; defaults to
    \dQuote{both}, but users can specify to show only the reads on
    \dQuote{plus} or \dQuote{minus} strand}
  \item{ylim}{range of read numbers to plot (y-axis limits); if not
    specified they are computed from the data in the specified region}
  \item{highlight}{currently unused}
  \item{main}{string; main title to use for the plot}
  \item{\dots}{further arguments passed on to the more primitive
    plotting functions used}
}
\details{
  This function implements the \code{plot} method for objects
  of class \code{AlignedGenomeIntervals}.
}
\value{
  Returns \code{NULL}; this function is called for the side-effect of
  creating the plot.
}
\author{Joern Toedling, Wolfgang Huber}
\note{
  This function was inspired by and borrows source code from the
  function \code{plotAlongChrom} in package \code{tilingArray}
}
\seealso{\code{\link{AlignedGenomeIntervals-class}}}
\examples{
## See the examples in the vignette and on the manual page
##  of the class 'AlignedGenomeIntervals'
}
\keyword{internal}
\keyword{hplot}
