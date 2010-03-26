\name{plotNegBinomFit}
\alias{plotNegBinomFit}
\title{Plot Negative Binomial Fit}
\description{
  Plot Negative Binomial Fit
}
\usage{
plotNegBinomFit(x, breaks = c(-0.5:15.5, 1e+08), bar.col=rainbow(2),
                addLegend=TRUE, legend.names=c("data","background"), ...)
}
\arguments{
  \item{x}{data.frame; slidingWindowSummary }
  \item{breaks}{numeric vector of breakpoints}
  \item{bar.col}{colours for the bars}
  \item{addLegend}{logical; should a legend be added in the top-right
    corner of the plot}
  \item{legend.names}{character; names for the legend}
  \item{\dots}{further arguments passed on to function barplot}
}
%\details{ }
\value{
  returns NULL; only called for its side-effect of producing the plot
}
%\references{}
\author{J Toedling}
\seealso{\code{\link{barplot}}}
%\examples{}
\keyword{hplot}
\keyword{internal}
