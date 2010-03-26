\name{plotReads}
\Rdversion{1.1}
\alias{plotReads}
\title{Function to plot aligned reads along the chromosome}
\description{
  Function to plot aligned reads along the chromosome
}
\usage{
plotReads(dat, ylim, strand = "plus", vpr, sampleColor = NULL,
   zeroLine = FALSE, main, pointSize = unit(1, "mm"),
   cexAxisLabel = 1, cexAxis = 1, ylab, ...)
}
\arguments{
  \item{dat}{a list with arguments
    \describe{
      \item{x.start}{integer; the genome start coordinates of the data
	to visualise}
      \item{x.end}{integer; the genome end coordinates of the data to
	visualise}
      \item{y}{numeric; the levels of the data to visualise}
      \item{flag}{numeric; specifies the category of each value,
	e.g. marks which data values belong to unique read alignments
	and which not}}}
  \item{ylim}{y-axis limits of the plotting window}
  \item{strand}{string; which of the two strands is plotted}
  \item{vpr}{which viewport to use for this plot}
  \item{sampleColor}{which colour to use for the data}
  \item{zeroLine}{logical; should a line at y=0 be drawn?}
  \item{main}{string; main title for the plot}
  \item{pointSize}{width of each dot/bar}
  \item{cexAxisLabel}{numeric; expansion factor for the axis labels}
  \item{cexAxis}{numeric; expansion factor for the axis labels}
  \item{ylab}{y-axis label}
  \item{\dots}{further arguments passed on to the more primitive
    plotting functions that are used}
}
\details{
  This function is used inside the plotting method for objects
  of class \code{AlignedGenomeIntervals}.
}
\value{
  returns \code{Null}; called for plotting single reads
  inside the function \code{plotAligned}
}
\author{Joern Toedling}
\keyword{hplot}
\keyword{internal}
