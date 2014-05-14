## Question: how to assess number of reads per window (function perWindow
##            and other use cases)

## Idea (already mentioned in cisGenome paper and analogous to estimation
##   of Poisson distribution parameters from limited count numbers in
##   association studies):
##  Negative binomial distribution as family of null distributions of
##  frequencies in uninteresting windows (cite CisGenome, BayesPeak)
##  Use observed count frequencies n0, n1, n2 of windows with 0, 1,
##  or 2 aligned reads for estimating the parameters of the negative
##  binomial distribution.
estimateNBParams <- function(n0, n1, n2){
  # build the ratios of these numbers to estimate the parameters:
  q1 <- n1/n0
  q2 <- n2/n1
  # then the estimates of the parameters are as follows:
  # a) dispersion parameter 'size' ('r' in Wikipedia notation)
  size.est <- q1/(2*q2-q1) # dispersion parameter size (or 'r'
  ## the mean 'mu' (or 'lambda' in Wikipedia notation),
  #mu.est <- q1*r.est/(r.est-q1)
  # which is equivalent to:
  mu.est <- q1/(1-2*q2+q1)
  return(list(size=as.vector(size.est), mu=as.vector(mu.est)))
}#estimateNBParams

## in some cases this Neg.bin estimation does not work,
### then use simpler Poisson
estimatePoissonParam <- function(n0, n1, n2){
  # build the ratios of these numbers to estimate the parameters:
  q1 <- n1/n0
  q2 <- n2/n1
  # then the estimates of the parameters are as follows:
  mu.est <- mean(q1, 2*q2)
  return(mu.est)
}#estimatePoissonParam

addNBSignificance <- function(x, estimate="NB.012", correct="none", max.n=10L){
  stopifnot(inherits(x, "slidingWindowSummary"),
            "n.reads" %in% names(x))
  correct  <- match.arg(correct, p.adjust.methods)
  estimate <- match.arg(estimate, c("NB.012", "NB.ML", "Poisson", "Poisson2"))
  ## 'x$n.reads' has to be integer, but maybe is not when weigthing each
  ##  interval count by the number of matches.
  x$n.reads <- round(x$n.reads)
  ## 1st way, using counts of windows with 0, 1, or 2 reads
  if (estimate == "NB.012"){
    tab.n <- table(x$n.reads)
    if (any(is.na(c(tab.n["0"], tab.n["1"], tab.n["2"]))))
      stop("Cannot estimate paramters of Negative-Binomial distribution, because there are no windows with 0, 1, or 2 aligned reads.")
    E <- estimateNBParams(tab.n["0"], tab.n["1"], tab.n["2"])
  }
  if (estimate == "NB.ML"){
    require("MASS")
    F <- suppressWarnings(MASS::fitdistr(x$n.reads[x$n.reads <= max.n],
                                         "negative binomial"))
    E <- as.list(F$estimate)
  }
  if (estimate == "Poisson"){
    tab.n <- table(x$n.reads)
    E <- list(mu=estimatePoissonParam(tab.n["0"], tab.n["1"], tab.n["2"]),
              size=.Machine$integer.max)
  }
  if (estimate == "Poisson2"){
    ## naive Poisson
    E <- list(mu=sum(x$n.reads)/nrow(x),
              size=.Machine$integer.max)
  }
  pval <- suppressWarnings(pnbinom(x$n.reads-1, mu=E$mu, size=E$size,
                                   lower.tail=FALSE))
  if (any(is.nan(pval)))
    stop("Parameter estimation '", estimate,"' failed. Try another method.\n")
  pval <- p.adjust(pval, method=correct)
  x$p.value <- pval
  E$"p.adjust.method" <- correct
  E$type <- estimate
  attr(x, "NBparams") <- E
  return(x)
}#addNBSignificance


#### function for plotting the fit:
plotNegBinomFit <- function(x, breaks=c(-0.5:15.5, 1e8),
                            bar.col=rainbow(2),
                            addLegend=TRUE,
                            legend.names=c("data","background"), ...)
{
  stopifnot(inherits(x, "slidingWindowSummary"),
            length(bar.col)==2,
            !is.null(attr(x, "NBparams")),
            all(c("size","mu") %in% names(attr(x, "NBparams"))))
  
  ## draw random read numbers (as many as in data)
  x.rand <- rnbinom(nrow(x), size=attr(x,"NBparams")$size,
                    mu=attr(x,"NBparams")$mu)
  xcut <- table(cut(x$n.reads, breaks=breaks))
  randcut <- table(cut(x.rand, breaks=breaks))
  mat <- rbind(xcut, randcut)
  barnames <- floor(breaks[-c(1,length(breaks))])
  barnames <- c(barnames, paste(">",barnames[length(barnames)],sep=""))
  barplot(mat, beside=TRUE, names.arg=barnames,
          ylab="Frequency", col=bar.col,...)
  if (addLegend)
    legend(x="topright", fill=bar.col, legend=legend.names)
  invisible(NULL)
}#plotNegBinomFit

# test:
#sel.cols <- brewer.pal(8, "Set2")[c(1,8)]
#plotNegBinomFit(resWins[[1]],xlab="Reads per window", bar.col=sel.cols)
