plotAllChrom <- function(x, ...){
  spacex <- chromosome(x)
  ## normalise chrosome names:
  spacex <- gsub("MT","M", spacex)
  spacex <- paste("chr",gsub("^chr","",spacex),sep="")
  allChr <- unique(spacex)
  nChr <- length(allChr)
  chrlens <- girafe:::getChromLengths(x)
  ## order chromosomes such that the combined length of two
  ##  chromosomes in each row is similar:
  if (nChr>1){
    par(mfrow=c(ceiling(nChr/2), 2))
    ## order chromosomes such that the combined length of two
    ##  chromosomes in each row is similar:
    ord <- order(chrlens, decreasing=TRUE)
    matord <- matrix(0, nrow=ceiling(nChr/2), ncol=2)
    matord[,1] <- ord[1:ceiling(nChr/2)]
    matord[1:floor(nChr/2),2] <- rev(ord[(ceiling(nChr/2)+1):nChr])
    ## matord gives order of chromosomes for layout()
    ##  now adjust width of each panel for layout()
    matChrlens <- matrix(0, nrow=ceiling(nChr/2), ncol=2)
    for (i in which(matord>0))
      matChrlens[i] <- chrlens[matord[i]]
    maxlen <- max(rowSums(matChrlens))
    matwidths <- t( apply(matChrlens, 1, function(z)
                          round(z*10 /sum(z))) )
    # t() because the result is transformed with apply
    #layout(matord, widths=matwidths, heights=rep(5, nrow(matord)))
  }
  for (thisChr in allChr){
    onChrom <- which(spacex==thisChr)
    these.x <- rowMeans(x[onChrom,1:2])
    these.y <- ifelse(strand(x)[onChrom]=="-",
                      -1*reads(x)[onChrom],
                      reads(x)[onChrom])
    plot(x=these.x, y=these.y, main=thisChr, type="h",
         xlab="Postion [bp]", ylab="Number of reads",
         xlim=c(1, chrlens[thisChr]),
         ylim=c(min(these.y, 0), max(these.y, 0)),
         frame.plot=FALSE, xaxt="n", ...)
    abline(h=0, lwd=2)
    these.ticks <- axTicks(side=1)
    nticks <- length(these.ticks)
    points(x=these.ticks, y=numeric(nticks), pch="|", cex=.7)
    text(x=these.ticks, y=numeric(nticks), labels=these.ticks,
         pos=ifelse(mean(these.y>0)>0.5, 1, 3), xpd=TRUE,
         font=2, cex=.7)
  }
  invisible(NULL)
} # plotAllChrom
