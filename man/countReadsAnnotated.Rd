\name{countReadsAnnotated}
\alias{countReadsAnnotated}
\title{Sum up aligned reads per category of genome feature}
\description{
  A function to sum up aligned reads per category of genome feature
  (i.e. gene, ncRNA, etc.).
}
\usage{
countReadsAnnotated(GI, M, typeColumn="type", fractionGI=0.7,
                    mem.friendly=FALSE, showAllTypes=FALSE)
}
\arguments{
  \item{GI}{object of class \code{AlignedGenomeIntervals}}
  \item{M}{Annotation object of class \code{Genome_intervals_stranded}
    or \code{Genome_intervals}; describes the genomic coordinates of
    annotated genome features, such as genes, miRNAs, etc.}
  \item{typeColumn}{string; which column of the annotation object
    \code{M} describes the type of the genome feature}
  \item{fractionGI}{which fraction of the intervals in object \code{GI}
    are required to ovelap with a feature in \code{M} in order to be
    considered to correspond to that feature.}
  \item{mem.friendly}{logical; should a version which requires less
    memory but takes a bit longer be used}
  \item{showAllTypes}{logical; should a table of genome feature types in
    \code{M} be displayed?}
}
\details{
  The read counts are summed up over each type of genome feature, and
  the read counts are normalised by their number of genomic matches. For
  example if a read has two matches in the genome, but only one inside a
  miRNA, it would count 0.5 for miRNAs.
}
\value{
  A named numeric vector which gives the summed read counts for each
  supplied type of genome feature.
}
\author{J Toedling}
%\seealso{\code{\link{help}}}
\examples{
  A <- AlignedGenomeIntervals(
           start=c(1,8,14,20), end=c(5,15,19,25),
           chromosome=rep("chr1", each=4),
           strand=c("+","+","+","+"),
           sequence=c("ACATT","TATCGGAC","TCGGACT","GTAACG"),
           reads=c(7L, 2L, 4L, 5L) )
  M2 <- new("Genome_intervals_stranded",
            rbind(c(2,6), c(1,15), c(20,30)),
            closed = matrix(TRUE, ncol=2, nrow=3),
            annotation = data.frame(
              seq_name= factor(rep("chr1", 3)),
              inter_base= logical(3),
              strand=factor(rep("+", 3), levels=c("+","-")),
              alias=c("miRNA1","gene1","tRNA1"),
              type=c("miRNA","gene","tRNA")) )
  if (interactive()){
     grid.newpage()
     plot(A, M2, chr="chr1", start=0, end=35,
          nameColum="alias", show="plus")
  }
  countReadsAnnotated(A, M2, typeColumn="type")
}
\keyword{manip}
