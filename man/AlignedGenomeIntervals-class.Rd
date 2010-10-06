\name{AlignedGenomeIntervals-class}
\Rdversion{1.1}
\docType{class}
\alias{AlignedGenomeIntervals-class}
\alias{AlignedGenomeIntervals}
\alias{[,AlignedGenomeIntervals,ANY,ANY-method}
\alias{c.AlignedGenomeIntervals}
\alias{c,AlignedGenomeIntervals-method}
\alias{clusters,AlignedGenomeIntervals-method}
\alias{clusters,Genome_intervals-method}
\alias{coerce,AlignedRead,AlignedGenomeIntervals-method}
\alias{coerce,AlignedGenomeIntervals,RangedData-method}
\alias{coverage,AlignedGenomeIntervals-method}
\alias{detail,AlignedGenomeIntervals-method}
\alias{extend}
\alias{extend,AlignedGenomeIntervals-method}
\alias{extend,Genome_intervals_stranded-method}
\alias{extend,Genome_intervals-method}
\alias{id,AlignedGenomeIntervals-method}
\alias{id<-}
\alias{id<-,AlignedGenomeIntervals,character-method}
\alias{hist,AlignedGenomeIntervals-method}
\alias{interval_included,AlignedGenomeIntervals,Genome_intervals_stranded-method}
\alias{interval_included,Genome_intervals_stranded,AlignedGenomeIntervals-method}
\alias{interval_included,AlignedGenomeIntervals,AlignedGenomeIntervals-method}
\alias{interval_overlap,AlignedGenomeIntervals,Genome_intervals-method}
\alias{interval_overlap,Genome_intervals,AlignedGenomeIntervals-method}
\alias{interval_overlap,AlignedGenomeIntervals,Genome_intervals_stranded-method}
\alias{interval_overlap,Genome_intervals_stranded,AlignedGenomeIntervals-method}
\alias{interval_overlap,AlignedGenomeIntervals,AlignedGenomeIntervals-method}
\alias{matches}
\alias{matches,AlignedGenomeIntervals-method}
\alias{matches<-}
\alias{matches<-,AlignedGenomeIntervals,integer-method}
\alias{nchar,AlignedGenomeIntervals-method}
\alias{organism,AlignedGenomeIntervals-method}
\alias{organism<-}
\alias{organism<-,AlignedGenomeIntervals,character-method}
\alias{plot,AlignedGenomeIntervals-method}
\alias{plot,AlignedGenomeIntervals,ANY-method}
\alias{plot,AlignedGenomeIntervals,missing-method}
\alias{plot,AlignedGenomeIntervals,Genome_intervals_stranded-method}
\alias{reads}
\alias{reads,AlignedGenomeIntervals-method}
\alias{reads<-}
\alias{reads<-,AlignedGenomeIntervals,character-method}
\alias{reduce,AlignedGenomeIntervals-method}
\alias{reduce,Genome_intervals-method}
\alias{sample,AlignedGenomeIntervals-method}
\alias{score,AlignedGenomeIntervals-method}
\alias{score<-,AlignedGenomeIntervals,numeric-method}
\alias{score<-}
\alias{seq_name}
\alias{seq_name,AlignedGenomeIntervals-method}
\alias{seq_name,Genome_intervals-method}
\alias{show,AlignedGenomeIntervals-method}
\alias{sort,AlignedGenomeIntervals-method}
\alias{strand,AlignedGenomeIntervals-method}
\alias{strand<-,AlignedGenomeIntervals,vector-method}
\alias{strand<-,AlignedGenomeIntervals,factor-method}
\alias{subset,AlignedGenomeIntervals-method}
\alias{width,AlignedGenomeIntervals-method}
\alias{export}
\alias{export,AlignedGenomeIntervals,character,character-method}
\alias{chromosome,AlignedGenomeIntervals-method}
\alias{chromosome,Genome_intervals-method}
\title{Class 'AlignedGenomeIntervals'}
\description{
  A class for representing reads from next-generation sequencing
  experiments that have been aligned to genomic intervals.
}
\section{Objects from the Class}{
  Objects can be created either by:
  \enumerate{
    \item calls of the form
    \code{new("AlignedGenomeIntervals", .Data, closed, ...)}.
    \item using the auxiliary function \code{AlignedGenomeIntervals} and
    supplying separate vectors of same length which hold the
    required information:\cr
    \code{AlignedGenomeIntervals(start, end, chromosome, strand, reads,
      matches, sequence)}\cr
    If arguments \code{reads} or \code{matches} are not specified, they
    are assumed to be '1' for all intervals.
    \item or, probably the most common way, by coercing from objects of
    class \code{AlignedRead}.
  }
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{two-column integer matrix, holding the
      start and end coordinates of the intervals on the chromosomes}
    \item{\code{sequence}:}{character; sequence of the read aligned to
      the interval}
    \item{\code{reads}:}{integer; total number of reads that were
      aligned to this interval}
    \item{\code{matches}:}{integer; the total number of genomic
      intervals that reads which were aligned to this interval were
      aligned to. A value of '1' thus means that this read sequence
      matches uniquely to this one genome interval only}
    \item{\code{organism}:}{string; an identifier for the genome of
      which organism the intervals are related to. Functions making use
      of this slot require a specific annotation package
      \code{org.<organism>.eg.db}. For example if \code{organism} is
      'Hs', the annotation package 'org.Hs.eg.db' is utilised by these
      functions. The annotation packages can be obtained from the
      Bioconductor repositories.}
    \item{\code{annotation}:}{data.frame; see class
      \code{genome_intervals} for details}
    \item{\code{closed}:}{matrix; see class
      \code{genome_intervals} for details}
    \item{\code{type}:}{character; see class
      \code{genome_intervals} for details}
    \item{\code{score}:}{numeric; optional score for each aligned genome
      interval}
    \item{\code{id}:}{character; optional identifier for each aligned
      genome interval}
  }
}
\section{Extends}{
Class \code{\link[genomeIntervals]{Genome_intervals-class}}, directly.
Class \code{\link[intervals:Intervals-class]{Intervals_full}}, by class
"Genome_intervals", distance 2.
}
\section{Methods}{
  \describe{
    \item{coerce}{Coercion method from objects of class
      \code{AlignedRead}, which is defined in package \code{ShortRead},
      to objects of class \code{AlignedGenomeIntervals}}
    \item{coerce}{Coercion method from objects of class
      \code{AlignedGenomeIntervals} to objects of class
      \code{RangedData}, which is defined in package \code{IRanges}}
    \item{coverage}{\code{signature("AlignedGenomeIntervals")}: computes
      the read coverage over all chromosomes. If the \code{organism} of
      the object is set correctly, the chromosome lengths are retrieved
      from the appropriate annotation package, otherwise the maximum
      interval end is taken to be the absolute length of that chromosome
      (strand).\cr
      The result of this method is a list and the individual list
      elements are of class \code{Rle}, a class for encoding long
      repetitive vectors that is defined in package \code{IRanges}.\cr
      The additional argument \code{byStrand} governs whether
      the coverage is computed separately for each strand. If
      \code{byStrand=FALSE} (default) only one result is returned per
      chromosome. If \code{byStrand=TRUE}, each chromosome result is
      again a list with two separate \code{Rle} objects.\cr
      By now, the \code{coverage} method for
      \code{AlignedGenomeIntervals} makes use of the method for
      \code{RangesList} objects from package \code{IRanges} (thanks to a
      suggestion from P. Aboyoun).
    }
    \item{detail}{\code{signature("AlignedGenomeIntervals")}: a more
      detailed output of all the intervals than provided by \code{show};
      only advisable for objects containing few intervals}
    \item{extend}{\code{signature("AlignedGenomeIntervals")} with
      additional arguments \code{fiveprime=0L} and
      \code{threeprime=0L}. These must be integer numbers and greater
      than or equal to 0. They specify how much is subtracted from the
      left border of the interval and added to the right side. Which end
      is 5' and which one is 3' are determined from the strand
      information of the object.
      Lastly, if the object has an \code{organism} annotation, it is
      checked that the right ends of the intervals do not exceed the
      respective chromosome lengths.}
    \item{export}{export the aligned intervals as tab-delimited text
      files which can be uploaded to the UCSC genome
      browser as \sQuote{custom tracks}.
      Currently, there are methods for exporting the data
      into \sQuote{bed} format and \sQuote{bedGraph} format,
      either writing the intervals from both strands into one file or
      into two separate files (formats \sQuote{bedStrand} and
      \sQuote{bedGraphStrand}, respectively).
      Details about these track formats can be found
      at the UCSC genome browser web pages.}
    \item{hist}{\code{signature("AlignedGenomeIntervals")}: creates
      a histogram of the lengths of the reads aligned to the intervals}
    \item{organism}{Get or set the organism that the genome intervals in
      the object correspond to. Should be a predefined code, such as
      'Mm' for mouse and 'Hs' for human. The reason for this code, that,
      if the organism is set, a corresponding annotation package that is
      called \code{org.<organism>.eg.db} is used, for example for
      obtaining the chromosome lengths to be used in methods such as
      \code{coverage}. These annotation packages can be obtained from
      the Bioconductor repository.
    }
    \item{plot}{visualisation method; a second argument of class
      \code{Genome_intervals_stranded} can be provided for additional
      annotation to the plot. Please see below and in the vignette for
      examples. Refer to the documentation of \code{\link{plotAligned}}
      for more details on the plotting function.}
    \item{reduce}{collapse/reduce aligned genome intervals by combining
      intervals which are completely included in each other, combining
      overlapping intervals AND combining immediately adjacent
      intervals. Intervals are only combined if they are on the same
      chromosome, the same strand AND have the same specificity of the
      aligned reads. \cr
      If you only want to combine intervals that have exactly the same
      start and stop position (but may have reads of slightly different
      sequence aligned to them), then use the optional argument
      \code{exact=TRUE}. \cr
      Finally, it's possible to only collapse/reduce aligned genome
      intervals that overlap each other by at least a certain fraction
      using the argument \code{min.frac}. \code{min.frac} is a number
      between 0.0 and 1.0. For example, if you call \code{reduce} with
      argument \code{min.frac=0.4}, only intervals that overlap
      each other by at least 40 percent are collapsed/merged.
    }
    \item{sample}{draw a random sample of \code{n} (Argument
      \code{size}) of the aligned reads (without or with replacement)
      and returns the \code{AlignedGenomeIntervals} object defined by
      these aligned reads.}
    \item{score}{access or set a custom score for the object}
    \item{sort}{sorts the intervals by chromosome name, start and end
      coordinate in increasing order (unless \code{decreasing=TRUE} is
      specified) and returns the sorted object}
    \item{subset}{take a subset of reads, matrix-like subsetting via
      '\[' can also be used}
    }
}
\author{Joern Toedling}
\seealso{
  \code{\link[genomeIntervals]{Genome_intervals-class}},
  \code{\link[ShortRead]{AlignedRead-class}},
  \code{\link[IRanges]{RangedData-class}},
  \code{\link[IRanges]{RangedData-class}},
  \code{\link{plotAligned}}
}
\examples{
  ############# toy example:
  A <- new("AlignedGenomeIntervals",
         .Data=cbind(c(1,3,4,5,8,10), c(5,5,6,8,9,11)),
         annotation=data.frame(
           seq_name=factor(rep(c("chr1","chr2","chr3"), each=2)),
           strand=factor(c("-","-","+","+","+","+") ,levels=c("-","+")),
           inter_base=rep(FALSE, 6)),
         reads=rep(3L, 6), matches=rep(1L,6),
         sequence=c("ACATT","ACA","CGT","GTAA","AG","CT"))

  show(A)
  detail(A)

  ## alternative initiation of this object:
  A <- AlignedGenomeIntervals(
     start=c(1,3,4,5,8,10), end=c(5,5,6,8,9,11),
     chromosome=rep(c("chr2","chrX","chr1"), each=2),
     strand=c("-","-","+","+","+","+"),
     sequence=c("ACATT","ACA","CGT","GTAA","AG","CT"),
     reads=rep(3L, 6))
  detail(A)

  ## custom identifiers can be assigned to the intervals
  id(A) <- paste("gi", 1:6, sep="")

  ## subsetting and combining
  detail(A[c(1:4)])
  detail(c(A[1], A[4]))

  ## sorting: always useful
  A <- sort(A)
  detail(A)
  
  ## the 'reduce' method provides a cleaned-up, compact set
  detail(reduce(A))
  ##  with arguments specifying additional conditions for merging
  detail(reduce(A, min.frac=0.8))

  ## 'sample' to draw a sample subset of reads and their intervals
  detail(sample(A, 10))
  
  ## biological example
  exDir <- system.file("extdata", package="girafe")
  exA   <- readAligned(dirPath=exDir, type="Bowtie", 
    pattern="aravinSRNA_23_no_adapter_excerpt_mm9_unmasked.bwtmap")
  exAI <- as(exA, "AlignedGenomeIntervals")
  organism(exAI) <- "Mm"
  show(exAI)
  ## which chromosomes are the intervals on?
  table(chromosome(exAI))

  ## subset
  exAI[is.element(chromosome(exAI),  c("chr1","chr2"))]

  ## compute coverage per chromosome:
  coverage(exAI[is.element(chromosome(exAI),  c("chr1","chr2"))])

  ### plotting:
  load(file.path(exDir, "anno_mm_genint.RData"))
  plot(exAI, mm.gi, chr="chrX", start=50400000, end=50410000)

  ### overlap with annotated genome elements:
  exOv <- interval_overlap(exAI, mm.gi)
  ## how many elements do read match positions generally overlap:
  table(listLen(exOv))
  ## what are the 12 elements overlapped by a single match position:
  mm.gi[exOv[[which(listLen(exOv)==12)]]]
  ## what kinds of elements are overlapped
  (tabOv <- table(as.character(mm.gi$type)[unlist(exOv)]))
  ### display those classes:
  my.cols <- rainbow(length(tabOv))
  pie(tabOv, col=my.cols, radius=0.95)
}
\keyword{classes}
