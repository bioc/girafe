perWindow <-
function(object, chr, winsize, step, normaliseByMatches=TRUE,
         mem.friendly=FALSE)
{
  ### arguments:
  # object: Modified genome_intervals object
  stopifnot(inherits(object, "AlignedGenomeIntervals"),
            is.logical(normaliseByMatches),
            is.character(chr), length(chr)==1L)
  object <- object[annotation(object)$"seq_name"==chr, ]
  if (nrow(object)==0){
    warning(paste("No aligned intervals on chromosome '", chr,"'",sep=""))
    return(NULL)
  }
  ## get vector of chromosome lengths
  chrLengths <- getChromLengths(object)
  stopifnot(chr %in% names(chrLengths))
  chrlen <- chrLengths[[chr]]

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
  ov <- interval_overlap(iv.win, object, mem.friendly=mem.friendly)

  ## read counts normalised by number of read matches or not?
  if (normaliseByMatches)
    sumfun <- function(x) sum(object@reads[x]/object@matches[x])
  else
    sumfun <- function(x) sum(object@reads[x])

  ## create results data frame
  res <- data.frame(chr=rep(chr, length(winstarts)),
                    start=winstarts, end=winends,
                    n.overlap=listLen(ov),
                    n.reads= sapply(ov, sumfun),
                    n.unique=sapply(ov, function(x)
                      sum(object@matches[x]==1)),
                    frac.plus=sapply(ov, function(x)
                      mean(strand(object[x])=="+")),
                    max.reads=sapply(ov, function(x){
                      xreads <- object@reads[x]
                      if (length(xreads)==0)
                        return(0)
                      else
                        return(max(object@reads[x]))}),
                    first=sapply(ov, function(x){
                      if (length(x)==0) return(as.integer(NA))
                      min(object[x, 1]) } ),
                    last=sapply(ov, function(x){
                      if (length(x)==0) return(as.integer(NA))
                      max(object[x, 2]) } )
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
                    # * 1/(last-first+1L))
  class(res) <- c(class(res), "slidingWindowSummary")
  return(res)
}#perWindow
