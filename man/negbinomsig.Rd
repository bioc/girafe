\name{addNBSignificance}
\alias{addNBSignificance}
\alias{estimateNBParams}
\title{assess significance of sliding-window read counts}
\description{
  This function can be used to assess the significance of sliding-window
  read counts. The background distribution of read counts in windows
  is assumed to be a Negative-Binomial (NB) one.
  The two parameters of the NB distribution, mean \sQuote{mu} and
  dispersion \sQuote{size}, are estimated using any of the methods
  described below (see details).
  The estimated NB distribution is used to assign a \emph{p}-value to
  each window based on the number of aligned reads in the window.
  The \emph{p}-values can be corrected for multiple testing using any
  of the correction methods implemented for \code{p.adjust}.
}
\usage{
addNBSignificance(x, estimate="NB.012", correct = "none",  max.n=10L)
}
\arguments{
  \item{x}{
    A \code{data.frame} of class \code{slidingWindowSummary}, as
    returned by the function \code{\link{perWindow}}.}
  \item{estimate}{
    string; which method to use to estimate the parameters of the NB
    background distribution; see below for details}
  \item{correct}{
    string; which method to use for \emph{p}-value adjustment;
    can be any method that is implemented for \code{p.adjust} including
    \dQuote{none} if no correction is desired.}
   \item{max.n}{integer; only relevant if \code{estimate=="NB.ML"}; in
     that case specifies that windows with up to this number of aligned
     reads should be considered for estimating the background
     distribution.}
}
\details{
  The two parameters of the Negative-Binomial (NB) distribution are:
  mean \sQuote{\eqn{\lambda}{lambda}} (or \sQuote{mu}) and size
  \sQuote{\eqn{r}} (or \sQuote{size}).\cr

  The function knows a number of methods to estimate the parameters of
  the NB distribution.
  \describe{
    \item{\dQuote{NB.012}}{Solely the windows with only 0,
      1, or 2 aligned reads are used for estimating
      {\eqn{\lambda}{lambda}} and \sQuote{\eqn{r}}.
      From the probability mass function \eqn{g(k)=P(X=k)} of the NB
      distribution, it follows that the ratios
      \deqn{q_1 = \frac{g(1)}{g(0)} = \frac{\lambda \cdot r}{\lambda +
	  r}}{q_1 = g(1)/g(0) = lambda r/(lambda+r)}
      and
      \deqn{q_2= \frac{g(2)}{g(1)}= \frac{\lambda \cdot (r+1)}{2 \cdot (\lambda + r)}~.}{q_2 = g(2)/g(1) = lambda (r+1)/(2 (lambda+r)).}
      The observed numbers of windows with 0-2 aligned reads are used to
      estimate
      \deqn{\widehat{q_1}=\frac{n_1}{n_0}}{q_1 = n_1/n_0}
      and
      \deqn{\widehat{q_2}=\frac{n_2}{n_1}}{q_2 = n_2/n_1}
      and from these estimates, one can obtain estimates for
      \eqn{\widehat{\lambda}}{'lambda'} and
      \eqn{\widehat{r}}{'r'}.}
    
    \item{\dQuote{NB.ML}}{This estimation method uses the function
      \code{fitdistr} from package \sQuote{MASS}. Windows with up to
      \code{n.max} aligned reads are considered for this estimate.}
    
    \item{\dQuote{Poisson}}{This estimate also uses the windows the 0-2
      aligned reads, but uses these numbers to estimates the parameter
      {\eqn{\lambda}{lambda}} of a Poisson distribution. The parameter
      \sQuote{\eqn{r}} is set to a very large number, such that the
      estimated NB distribution actually is a Poisson distribution with
      mean and variance equal to \eqn{\lambda}{lambda}.}
  }
}
\value{
  A \code{data.frame} of class \code{slidingWindowSummary}, which is the
  the supplied argument \code{x} extended by an additional column
  \code{p.value} which holds the \emph{p}-value for each window.
  The attribute \code{NBparams} of the result contains the list of the 
  estimated parameters of the Negative-Binomial background
  distribution.
}
\references{
  Such an estimation of the Negative-Binomial parameters has also been
  described in the paper:\cr
  Ji et al.(2008) An integrated system CisGenome for analyzing
  ChIP-chip and ChIP-seq data. Nat Biotechnol. 26(11):1293-1300.
}
\author{Joern Toedling}
\seealso{
  \code{\link{perWindow}}, \code{\link[stats]{p.adjust}}
}
\examples{
 exDir <- system.file("extdata", package="girafe")
 exA   <- readAligned(dirPath=exDir, type="Bowtie",
   pattern="aravinSRNA_23_no_adapter_excerpt_mm9_unmasked.bwtmap")
 exAI  <- as(exA, "AlignedGenomeIntervals")
 exPX  <- perWindow(exAI, chr="chrX", winsize=1e5, step=0.5e5)
 exPX  <- addNBSignificance(exPX, correct="bonferroni")
 str(exPX)
 exPX[exPX$p.value <= 0.05,]
}
\keyword{manip}
