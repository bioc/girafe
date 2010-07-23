\name{getFeatureCounts}
\alias{getFeatureCounts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{get the read counts for a supplied set of genomic features}
\description{
  get the read counts for a supplied set of genomic features
}
\usage{
getFeatureCounts(AI, FG, nameColumn = "Name", fractionIncluded = 1,
returnType = "AlignedGenomeIntervals", normalise = TRUE)
}
\arguments{
  \item{AI}{AlignedGenomeIntervals object}
  \item{FG}{Genome_intervals objects of genomic features}
  \item{nameColumn}{
%%     ~~Describe \code{nameColumn} here~~
}
  \item{fractionIncluded}{
%%     ~~Describe \code{fractionIncluded} here~~
}
  \item{returnType}{one of \code{AlignedGenomeIntervals} or
    \code{integer}}
  \item{normalise}{logical}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{Joern Toedling}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
\keyword{internal}% __ONLY ONE__ keyword per line
