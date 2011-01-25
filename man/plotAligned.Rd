\name{plotAligned}
\Rdversion{1.1}
\alias{plotAligned}
\title{Visualise reads aligned to genome intervals}
\description{
  Visualise reads aligned to genome intervals
}
\usage{
plotAligned(x, y, chr, start, end, plus.col = "#00441b",
   minus.col = "#283d78", gff, featureLegend = FALSE,
   gffChrColumn = "seq_name", gffTypeColumn="type",
   gffNameColumn="ID",
   featureExclude = c("chromosome", "nucleotide_match", "insertion"),
   showStrands="both", extraColors=NULL, ylim, highlight, main, ...)
}
\arguments{
  \item{x}{Object of class \code{AlignedGenomeIntervals}}
  \item{y}{This argument is only specified for compatibility with
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
    used to further annotate the plot. Note that it must include a
    column called \dQuote{type} that indicates the type of each genomic
    feature (e.g. miRNA, gene etc.). }
  \item{featureLegend}{logical; should a legend that describes the
    colour code for the annotated genome features be appended at the
    bottom of the plot?}
  \item{gffChrColumn}{string; which column of the \code{gff} data.frame
    holds the chromosome identifier of each feature.}
  \item{gffTypeColumn}{string; which column of the \code{gff} data.frame
    holds the type/class identifier of each feature. Used for the
    colouring of features.}
  \item{gffNameColumn}{what is the column of the gff data.frame called
    that holds the identifier of the element that should be displayed in
    the plot; default: \dQuote{name}}
  \item{featureExclude}{character; which kinds of annotated genome
    features specified in the \code{gff} are to be ignored for the plot}
  \item{showStrands}{string; which strands to show in the plot; defaults
    to \dQuote{both}, but users can specify to show only the reads on
    \dQuote{plus} or \dQuote{minus} strand}
  \item{extraColors}{named character vector which allows the user to
    specify custom colours for feature types; colours must be specified
    in RGB format as hexadecimal strings starting with \dQuote{#},
    e.g. \dQuote{#addfff} for light-blue}
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
  A <- AlignedGenomeIntervals(
        start=c(1,8,10,20), end=c(5,15,16,25),
        chromosome=rep("chr1", each=4),
        strand=c("+","+","+","+"),
        sequence=c("ACATT","TATCGGAC","TCGGACT","GTAACG"),
        reads=c(5L, 2L, 4L, 7L) )
  M2 <- new("Genome_intervals_stranded",
            rbind(c(2,6), c(1,15), c(20,30)),
            closed = matrix(TRUE, ncol=2, nrow=3),
            annotation = data.frame(
              seq_name= factor(rep("chr1", 3)),
              inter_base= logical(3),
              strand=factor(rep("+", 3), levels=c("+","-")),
              alias=c("miRNA1","gene1","tRNA1"),
              type=c("miRNA","gene","tRNA")) )
  if (interactive())
    plot(A, M2, chr="chr1", start=0, end=35,
         nameColum="alias", showStrands="plus")

  if (interactive())
  ## use 'extraColors' to replace default colours (or to add new ones):
    plot(A, M2, chr="chr1", start=0, end=35, nameColum="alias",
         showStrands="plus", extraColors=c("tRNA"="#6666DD"))

  ## See also the examples in the vignette and on the manual page
  ##  of the class 'AlignedGenomeIntervals'
}
\keyword{internal}
\keyword{hplot}
