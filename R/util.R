###################################################################
### auxiliary functions:
###################################################################

# get a named integer vector of chromosome lengths, based on the
#  organism slot of AlignedGenomeIntervals objects
getChromLengths <- function(x){
  stopifnot(inherits(x, "AlignedGenomeIntervals"))
  if ( length(x@organism) != 1 )
    stop("Organism of '",deparse(substitute(x)),
         "' not defined! Set with organism <-")
  orgpackage <- paste("org",x@organism,"eg.db", sep=".")
  worked <- require(orgpackage, character.only=TRUE)
  if (!worked)
    stop("No package called '",orgpackage,"' found.",
         "Install this package or check wether the organism annotation",
         "of ",deparse(substitute(x))," is correct (e.g. 'Mm' or 'Hs')\n")
  chrlens <- get(paste(gsub("\\.db$", "", orgpackage),
                       "CHRLENGTHS", sep=""))
  names(chrlens) <- paste("chr",
                          gsub("^chr","", names(chrlens)),sep="")
  return(chrlens)
} #getChromLengths
