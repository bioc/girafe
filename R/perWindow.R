perWindow <-
function(object, chr, winsize, step, organism="Mm"){

  ### arguments:
  # object: Modified genome_intervals object
  stopifnot(inherits(object, "AlignedGenomeIntervals"))
  object <- object[annotation(object)$"seq_name"==chr, ]
  if (nrow(object)==0){
    warning(paste("No aligned intervals on chromosome '", chr,"'",sep=""))
    return(NULL)
  }
  ## organism package for chr length
  orgname <- paste("org.", organism, ".eg",sep="")
  require(paste(orgname, "db", sep="."), character.only=TRUE)
  stopifnot(exists(paste(orgname, "CHRLENGTHS",sep="")))
  chrlen <- get(paste(orgname, "CHRLENGTHS",sep=""))[gsub("^chr","",gsub("MT","M",chr))]
  stopifnot(!is.na(chrlen))

  ### prepare sliding windows
  winstarts <- seq.int(from=min(object[,1]), to=chrlen, by=step)
  winends <- winstarts + winsize - 1L
  if (winends[length(winends)]> chrlen){
    winends[length(winends)] <- chrlen
    winstarts[length(winends)] <- max(1, chrlen-winsize+1)
  } # if last window extends beyond chromosome end

  ### intervals for windows:
  iv.win <- new("Genome_intervals",
                cbind(winstarts, winends),
                annotation=data.frame(
                  "seq_name"=factor(rep(chr, length(winstarts))),
                  inter_base=vector("logical", length(winstarts)))
                )

  ### determine overlap:
  ov <- interval_overlap(iv.win, object)

  res <- data.frame(chr=rep(chr, length(winstarts)),
                    start=winstarts, end=winends,
                    n.overlap=listLen(ov),
                    n.reads= sapply(ov, function(x)
                      sum(object@reads[x])),
                    n.unique=sapply(ov, function(x)
                      sum(object@matches[x]==1)),
                    frac.plus=sapply(ov, function(x)
                      mean(strand(object[x])=="+")),
                    max.reads=sapply(ov, function(x){
                      xreads <- object@reads[x]
                      if (length(xreads)==0)
                        return(0)
                      else
                        return(max(object@reads[x]))})
                    )

  ### we want a high cluster score for
  ## a.) many match positions
  ## b.) few reads per match position
  ## c.) higher proportion of unique matches in the region
  ## d.) shorter region length
  ## e.) from both strands
  #res$score <- with(res, n.overlap*
  #                  (n.overlap/n.reads + n.unique/n.overlap -
  #                  (frac.plus-0.5)^2))
                    # * 1/(end-start+1))
  class(res) <- c(class(res), "slidingWindowSummary")
  return(res)
}#perWindow
