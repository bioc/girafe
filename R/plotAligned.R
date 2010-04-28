### visualization function for AlignedGenomeIntervals objects
plotAligned <- function(x, y, chr, start, end,
                        plus.col="#00441b", minus.col="#283d78",
                        gff, featureLegend=FALSE, gffChrColumn="seq_id",
                        gffNameColumn="name",
                        featureExclude=c("chromosome", "nucleotide_match",
                          "insertion"), show="both",
                        ylim, highlight, main,...){
  ### evaluate arguments:
  show <- match.arg(show, c("both", "plus", "minus"))
  xlim <- c(start, end)
  if (!missing(ylim))
    stopifnot(is.numeric(ylim), length(ylim)==2)
  thisCall <- match.call()
  if (!is.element(chr, levels(seq_name(x))))
    warning(paste("'",deparse(substitute(chr)),
                  "' is not among the chromosome/sequences names in ",
                  "this object, which are:\n",
                  paste(levels(seq_name(x)), collapse=","), sep=""))
  x <-  x[as.character(seq_name(x))==chr & (x[,1]>=start) & (x[,2]<=end)]
  if (nrow(x)==0)
    warning("No intervals with aligned reads in specified region.\n")
  x.plus <- x[annotation(x)$strand == "+"]
  x.minus <- x[annotation(x)$strand == "-"]
  max.plus <- suppressWarnings(max(1, max(x.plus@reads)))
  max.minus <- suppressWarnings(max(1, max(x.minus@reads)))
  
  ##-----------------------------------------------#
  ## PART II: plotting
  ##-----------------------------------------------#
  ## set up the viewports of the plot layout.
  VP <- c("title"=0.4, "readsplus"=5, "gff+"=1,
          "coord"=1, "gff-"=1, "readsminus"=5, "legend"=1)

  ## if desired show only one of the two strands
  if (show=="minus") VP <- VP[-match(c("readsplus", "gff+"), names(VP))]
  if (show=="plus")  VP <- VP[-match(c("readsminus", "gff-"), names(VP))]

  if(!featureLegend)
    VP <- VP[-which(names(VP)=="legend")]
  
  if(missing(gff))
    VP <- VP[-which(names(VP)=="gff+" | names(VP)=="gff-")]

  ## set default colours
  defaultColors <- c("plus"=plus.col, "minus"=minus.col,
                     "duplicated"="grey", "highlight"="red")
  
  ## plot margin
  pushViewport(viewport(width=0.85, height=0.95)) ## plot margin
  pushViewport(viewport(layout=grid.layout(length(VP), 1, height=VP)))

  names.strand <- c("plus"="Watson", "minus"="Crick")
  if (show!="both")
    show.strands <- show
  else
    show.strands <- c("plus", "minus")

  for (strand in show.strands){
    ylab <- paste("Reads on",names.strand[strand],"strand")
    x.strand <- get(paste("x", strand, sep="."))
    dat <- list(x.start = x.strand[,1],
                x.end =x.strand[,2],
                y = x.strand@reads,
                flag = as.numeric(x.strand@matches>1)*3)
    ## At this point, no matter what the input to the function was,
    ## we have the list 'dat' with elements
    ## x.start: x0 coordinate
    ## x.end : x1 coordinate
    ## y: y coordinate
    ## flag
              
    ## set up the y-axis limits
    if (missing(ylim))
      ylimdata <- c(0, get(paste("max",strand,sep=".")))
    else
      ylimdata <- ylim
      
    if (strand=="minus")
      ylimdata <- rev(-1*ylimdata)
    
    ## plot the data
    vpr <- which(names(VP)==paste("reads", strand, sep=""))
    
    ######################################################
    ### set up viewport for plotting aligned reads #######
    ######################################################
    pushViewport(dataViewport(xData=xlim, yData=ylimdata,
                              extension=0, clip="off",
                              layout.pos.col=1, layout.pos.row=vpr))
    ypos <- unique(round(pretty(ylimdata)))
    grid.yaxis(at=ypos, label=abs(ypos), gp=gpar(cex=0.8),...)
    grid.text(ylab, x=-0.075, y=0.5,#mean(ylimdata),
              just=c("center", "center"),#"center"),
              default.units="npc", rot=90,
              gp=gpar(cex=0.9, font=2),...)

    ## only if there are reads on these strand:
    if (nrow(x.strand)>0){
      pushViewport(dataViewport(xData=xlim, yData=ylimdata, extension=0,
                                clip="on", layout.pos.col=1,
                                layout.pos.row=vpr))
      plotReads(dat, sampleColor=defaultColors[strand],
                strand=strand, vpr=vpr, ylim=ylimdata,...)
      popViewport(1)
    }
    ### end plotting the reads of the strand
    popViewport(1)
  } # for strand in c("plus","minus")
  
  ########################################################
  ### plot annotated genome features (supplied in the gff) 
  ########################################################
  if(!missing(gff))
    for (strand in c("plus"="+","minus"="-")[show.strands])
      plotFeatures(gff=gff, chr=chr, xlim=c(start, end),
                   strand=strand,
                   featureExclude=featureExclude,
                   featureColorScheme=1,
                   vpr=which(names(VP)==sprintf("gff%s", strand)),
                   gffChrColumn=gffChrColumn,
                   gffNameColumn=gffNameColumn, ...)
  
  ###########################################################
  ### plot  chromosomal coordinate axis #####################
  ###########################################################
  pushViewport(dataViewport(xData=xlim, yscale=c(-0.4,0.8),
                            extension=0, layout.pos.col=1,
                            layout.pos.row=which(names(VP)=="coord")))
  grid.lines(xlim, c(0,0), default.units = "native")
  tck <- alongChromTicks(xlim)
  grid.text(label=formatC(tck, format="d"), x=tck, y=0.2,
            just=c("centre", "bottom"), gp = gpar(cex=.6),
            default.units="native")
  grid.segments(x0 = tck, x1 = tck, y0 = -0.17, y1 = 0.17,
                default.units = "native")
            
  #### highlight (not implemented properly yet)
  if(!missing(highlight)){
    mt = (match(highlight$strand, c("-", "+"))-1.5)*2
    if (!is.null(highlight$coord))
      co <- highlight$coord
    else
      co <- highlight$x
    if(is.na(mt) || !is.numeric(co))
      stop("Invalid parameter 'highlight'.")
    strand.num <- ifelse(highlight$strand=="-",-1,1)
    grid.segments(x0=co, x1=co+(500*strand.num),
                  y0=c(0.4,0.4)*mt, y1=c(0.4,0.4)*mt,
                  default.units = "native", arrow=arrow(),
                  gp=gpar(col="violetred4", lwd=4))
  }

  ## end coordinate axis
  popViewport()
            
  #### TITLE
  pushViewport(viewport(layout.pos.col=1,
                        layout.pos.row=which(names(VP)=="title")))
  grid.text(label=paste("Chromosome", chr, "coordinate [bp]"),
            x=0.5, y=1, just="centre", gp=gpar(cex=1, font=2))
  if(!missing(main))
    grid.text(label=main, x=0.05, y=1, just="centre",
              gp=gpar(cex=1, font=2))
  popViewport()
  
  #### FEATURE LEGEND
  if(featureLegend)
    plotAlongChromLegend(which(names(VP)=="legend"),
                         featureColorScheme=1,
                         featureExclude=featureExclude)

  #### END ALL PLOTTING
  popViewport(2)
  invisible(dat)
  
} # plotAligned


##### ------------------------------------------------------------
## auxiliary function: plot Features
## ------------------------------------------------------------
plotFeatures <- function(gff, chr, xlim, strand, vpr, featureColorScheme=1,
                         featureExclude=c("chromosome", "nucleotide_match",
                           "insertion"), featureNoLabel=c("uORF", "CDS"),
                         gffNameColumn="name",
                         gffChrColumn="seq_id", ...)
{
  #### check arguments:
  stopifnot(is.data.frame(gff),
            all(c(gffNameColumn, gffChrColumn,
                  "strand","start","end") %in% names(gff)),
            length(strand)==1, strand %in% c("+", "-"))
  
  translateStrand <- function(x){
    if (x==-1) return("-")
    if (x==1) return("+")
    return(NA)}
  stopifnot(all(gff[,"start"] <= gff[, "end"]))
  if (is.numeric(gff$strand))
    gff$strand <- sapply(gff$strand, translateStrand)
   
  pushViewport(dataViewport(xData=xlim, yscale=c(-1.2,1.2),  extension=0,
                            clip="on", layout.pos.col=1, layout.pos.row=vpr))

  sel <- which(gff[, gffChrColumn]     == chr &
               gff[, "strand"]  == strand &
               gff[, "start"]   <= xlim[2] &
               gff[, "end"]     >= xlim[1])

  ## for label, use "symbol" if available, otherwise "name"
  featName = gff[sel, gffNameColumn]
  
  ## split by feature type (e.g. CDS, ncRNA)
  feature  = as.character(gff[sel, "type"])
  featsp = split(seq(along=sel), feature)
  
  ## There are now five different cases, and we need to deal with them:
  ## - ignorable features, given by featureExclude
  ## - genes: a horizontal line + name
  ## - introns: a caret
  ## - CDS: a box + no name
  ## - all others: a colored box + name

  ## in this vector we save those features for which we want to have names
  whnames = integer(0)

  ## 1. drop the ignorable ones
  featsp = featsp[ ! (names(featsp) %in% featureExclude) ]
  
  ## 3.introns
  wh = ("intron" == names(featsp))
  if(any(wh)) {
    i = featsp[["intron"]]
    s = sel[i]
    mid = (gff$start[s]+gff$end[s])/2
    wid = (gff$end[s]-gff$start[s])/2 
    for(z in c(-1,1))
      grid.segments(x0 = mid,
                    x1 = mid+z*wid,
                    y0 = 1.20*c("+"=1, "-"=-1)[strand],  ## istrand is 1 or 2
                    y1 = 0.95*c("+"=1, "-"=-1)[strand],
                    default.units = "native",
                    gp = gpar(col="black"))
     featsp = featsp[!wh]
  } ## if
  
  ## 4. colors for boxes
  ## check that we know how deal with all features
  featCols = featureColors(featureColorScheme)
  whm = names(featsp) %in% rownames(featCols)
  ### indicate features of unknown types as such:
  if(!all(whm)){
    warning("Don't know how to handle feature of type(s) '",
            paste(names(featsp)[!whm], collapse=", "), "' in gff.", sep="")
    names(featsp)[!whm] <- "unknown"
  }

  sfeatsp  = featsp[rownames(featCols)]
  ll       = listLen(sfeatsp)
  
  if(any(ll>0)) {
    i  = unlist(sfeatsp)
    gp = gpar(col = rep(featCols$col,  ll),
                 fill = rep(featCols$fill, ll))
    s  = sel[i]
    grid.rect(x     = gff$start[s],
              y     = 0,
              width = gff$end[s]-gff$start[s],
              height= 2,
              default.units = "native",
              just  = c("left", "center"),
              gp    = gp)
    whnames = c(whnames, unlist(sfeatsp[!(names(sfeatsp) %in% featureNoLabel)]))
    ## additional potentially useful values for featureNoLabel: "binding_site", "TF_binding_site"
  }

  ## labels
  if( !all(tolower(featureNoLabel)=="all") && (length(whnames)>0)) {

    ## this is a bit of a hack to abbreviate the labels of
    ##  "binding site" features:
    bindingRegexpr = "binding.?site.*$"
    isBindingSite = (regexpr(bindingRegexpr, featName[whnames]) > 0)
    if(any(isBindingSite)) {
      ## replace long labels
      featName[whnames] = gsub(bindingRegexpr, "bs", featName[whnames])
    }

    ## remove duplicated names that are not binding sites
    whnames = whnames[isBindingSite | !duplicated(featName[whnames])]

    txtcex = 0.6
    txtdy  = 0.7
    s      = sel[whnames]
    txtx   = (pmax.int(min(xlim), gff$start[s])+ pmin.int(max(xlim), gff$end[s]))/2
    txty   = numeric(length(s))
    ord    = order(txtx)
    whnames = whnames[ord]
    s      = s[ord]
    txtx   = txtx[ord]
    
    strw   = convertWidth(stringWidth(featName[whnames]), "native", valueOnly=TRUE)*txtcex
    rightB = txtx[1] + 0.5*strw[1]
    doText = rep(TRUE, length(whnames))

    # adjust text labels to be still readable in feature-dense areas:
    if(length(whnames) >1) {
      for(k in 2:length(whnames)) {
        leftB = txtx[k] - 0.5*strw[k]
        if(leftB > rightB) { # all texts not overlapping next to each other?
          rightB = txtx[k] + 0.5*strw[k]
        } else { # any overlaps?
          if(!any(txty[k-(1:2)]==txtdy)) {#  2 previous labels not moved up?
            txty[k]= txtdy                #   then this one 
          } else {                        #  else try move down:
            if(!any(txty[k-(1:2)]== -txtdy)) { 
              txty[k]= -txtdy             #  if 2 previous ones weren't
            } else {
              doText[k] = FALSE           #  otherwise don't put the label
            }
          }
        } ##  else
      } ## for
    }
    grid.text(label = featName[whnames][doText],
              x = txtx[doText], y = txty[doText], gp=gpar(cex=txtcex), 
              default.units = "native")
  } ## if
  popViewport()
} ## plotFeatures

##------------------------------------------------------------
##
##------------------------------------------------------------
alongChromTicks <- function(x){
  rx = range(x)
  lz = log((rx[2]-rx[1])/3, 10)
  fl = floor(lz)
  if( lz-fl > log(5, 10))
    fl = fl +  log(5, 10)
  tw = round(10^fl)
  i0 = ceiling(rx[1]/tw)
  i1 = floor(rx[2]/tw)
  seq(i0, i1)*tw
}# alongChromTicks

##------------------------------------------------------------
## featureColors
## note that features are drawn in the order in which they appear
## here, this can be used to let important features overdraw less
## important ones (e.g. tRNA is more specific than ncRNA)
## to test, say tilingArray:::plotAlongChromLegend()
##------------------------------------------------------------
featureColors <- function(scheme=1){
  defaultColors = c(
    "centromere"               = "#FFEDA0",    ## orange
    "telomere"                 = "#FFEDA0",    ## orange
    "novel_pseudogene"         = "#f0f0f0",    ## light gray
    "novel_retrotransposed"    = "#f0f0f0",    ## light gray
    "novel_coding"             = "#e0f1f2",    ## lighter blue
    "novel_RNA"                = "#b6e97a",    ## green 
    "novel_tRNA"               = "#b6e97a",    ## green 
    "novel_scRNA"              = "#9C9BC1",    ## purple
    "novel_snRNA"              = "#9C7BC1",    ## purple
    "novel_rRNA"               = "#ffbe71",    ## meat
    "novel_snoRNA"             = "#8F6A68",    ## red brown
    "novel_miRNA"              = "#DC76DC",    ## light red-violet
    "pseudogene"               = "#e0e0e0",    ## light gray
    "uORF"                     = "#FED976" ,   ## orange
    "nc_primary_transcript"    = "#a0a0a0",    ## grey
    "repeat_family"            = "#CC6666",    ## light red
    "repeat_region"            = "#e31a1c",    ## bright red
    "retrotransposed"          = "#f1b6da",    ## pink
    "transposable_element"     = "#f1b6da",    ## pink
    "transposable_element_gene"= "#f1b6da",
    "ARS"         = "#CC9966",             ## light brown
    "insertion"   = "#FFEDA0",             ## orange
    "CDS_dubious" = "#e0f1f2",             ## lighter blue
    "gene"        = "#addfff",             ## light blue
    "CDS"         = "#addfff",             ## light blue
    "coding"      = "#5E88B0",             ## light blue
    "exon"        = "#5E88B0",             ## aquamarine
    "transcript"  = "#5E88B0",             ## aquamarine
    "ncRNA"       = "#86b94a",             ## dark green 
    "tRNA"        = "#a6d96a",             ## green
    "snRNA"       = "#8C6BB1",             ## purple
    "rRNA"        = "#fdae61",             ## meat
    "ribozyme"    = "#dd8e41",             ## meat
    "snoRNA"      = "#7F5A58",             ## red brown
    "miRNA"       = "#cc66cc",             ## light red-violet
    "unknown"     = "#a0a0a0",
    "nucleosome_binding_motif" = "#C9C299",## lemon chiffon
    "TF_binding_site" = "#C9C299"          ## lemon chiffon
    )
  darkenborder <- as.logical(c(rep(1, length(defaultColors)-2),0, 0))
  stopifnot(length(darkenborder)==length(defaultColors))
  
  fill = switch(scheme,
    default  = defaultColors,
    unicolor = ifelse(is.na(defaultColors), NA,  "#addfff"),  ## light blue
    stop("Encountered error when filling in colors."))
  
  ## calculate hex string for a color that is a little bit darker than the
  ## hex string in the argument
  darken <- function(x, factor=0.5) {
    wh = which(!is.na(x))
    hex = sapply(x[wh], substring, first=c(2,4,6), last=c(3,5,7))
    hex = apply(hex, 2, function(h) as.integer(factor*as.integer(paste("0x", h, sep=""))))
    res = rep(as.character(NA), length(x))
    res[wh] = apply(hex, 2, function(h) sprintf("#%02x%02x%02x", h[1], h[2], h[3]))
    return(res)
  }
  border <- ifelse(darkenborder, darken(fill), fill)
  res <- data.frame(fill=I(fill), col =I(border))
  rownames(res) <- names(defaultColors) 
  return(res)
} # function featureColors

##------------------------------------------------------------
## legend
##------------------------------------------------------------
plotAlongChromLegend <- function(vpr, nr=2, featureColorScheme=1,
    featureExclude=c("chromosome", "nucleotide_match", "insertion"),
    mainLegend, cexLegend=0.35, cexMain=1)
{
  endVP = FALSE
  # when this function is called on its own 
  # set up a viewport
  if(missing(vpr)) { 
     endVP=TRUE      
     vpr = newVP(main=mainLegend, dataPanelHeight=1, cexMain=cexMain) # newVP sets up a new viewport
  }
  formatRow = function(featColsOneRow, row) {
    ## print(featColsOneRow)
    strWid   = convertWidth(stringWidth(rownames(featColsOneRow)), "npc", valueOnly=TRUE)
    n        = length(strWid)
    inbetWid = 0.2*min(strWid)
    totWid   = sum(strWid)+(n-1)*inbetWid
    x        = c(0, cumsum(strWid[-n])) + (0:(n-1))*inbetWid 
    y        = numeric(length(x))

    x      = x/totWid
    strWid = strWid/totWid
    grid.rect(x = x, width = strWid, 
              y = unit(row, "native"), height = unit(1, "native")- unit(1, "mm"), 
              just  = c("left", "center"), default.units="npc",
              gp    = do.call(gpar, featColsOneRow))
    
    grid.text(label = rownames(featColsOneRow),
              x = unit(x + strWid/2, "native"), y = unit(row, "native"),
              just  = c("center", "center"), gp=gpar(cex=cexLegend))
  } 

  featCols = featureColors(featureColorScheme)
  featCols = featCols[ !(rownames(featCols) %in% featureExclude), ]

  pushViewport(viewport(layout.pos.col=1, layout.pos.row=vpr, yscale=c(0.5, nr+0.5)))

  i = 1:nrow(featCols)
  for(r in 1:nr)
    formatRow(featCols[ceiling(i/nrow(featCols)*nr-1e-10)==r, ], row=nr-r+1)
  
  popViewport()

  if(endVP)
     popViewport(2)
}


# this function sets up a new viewport.  It is used by plotAlongChromLegend, 
# plotSegmentationHeatmap and plotOneChIPSample when they are called as 
# stand-alone functions (ie when vpr is not specified)
newVP <- function(main, cexMain=1, dataPanelHeight=1, vpHeight=0.7, titleOffSet=0) {
  if(!missing(main)) {
    vpr = c("title"=0.1, "data"=dataPanelHeight)
    pushViewport(viewport(width=0.85, height=vpHeight)) ## plot margin
    pushViewport(viewport(layout=grid.layout(length(vpr), 1, height=vpr)))  
    pushViewport(viewport(layout.pos.col=1, layout.pos.row=which(names(vpr)=="title")))
    grid.text(label=main, x=0.5, y=1.1+titleOffSet, just="centre", gp=gpar(cex=cexMain))  
    popViewport()
    vpr = which(names(vpr)=="data")
  } else {
    vpr = c("data"=dataPanelHeight)
    pushViewport(viewport(width=0.85, height=vpHeight)) ## plot margin
    pushViewport(viewport(layout=grid.layout(length(vpr), 1, height=vpr)))
  }
  return(vpr)
}#newVP



##-------------------------------------------------------------
##   plot one ChIP-chip sample with lines and dots
##-------------------------------------------------------------
plotReads <- function(dat, ylim, strand="plus", vpr, sampleColor=NULL,
                      zeroLine=FALSE, main, pointSize=unit(1.0, "mm"),
                      cexAxisLabel=1, cexAxis=1, ylab, ...)
{
  endVP = FALSE
  ## in case this function is not called from plotAlignedGenomeIntervals
  if(missing(vpr)) {
     endVP=TRUE
     vpr = newVP(main=main, dataPanelHeight=1, vpHeight=0.95,
       titleOffSet=-0.9)
  }

  stopifnot(length(dat$y)==length(dat$x.start),
            length(dat$flag)==length(dat$x.start),
            length(dat$x.end)==length(dat$x.start))
  xorg  = dat$x
  
  if(missing(ylim)){
    suppressWarnings(ylim <- c(0, max(1, max(dat$y))))
    if (strand == "minus") ylim <- rev(-1*ylim)
  }
  xlimdata <- suppressWarnings(range(c(dat$x.start, dat$x.end)))

  ## the reads data. use two viewports for different clipping behavior
  if (missing(vpr)){
    if(!missing(ylab)) {
      pushViewport(dataViewport(xData=xlimdata, yData=ylim, extension=0, clip="off",
                                layout.pos.col=1, layout.pos.row=vpr))
      grid.yaxis(gp=gpar(cex=cexAxis),...)
      grid.text(ylab, x=-0.075, y=0.5, rot=90, gp=gpar(cex=cexAxisLabel),...)
      pushViewport(dataViewport(xData=xlimdata, yData=ylim, extension=0, clip="on",
                                layout.pos.col=1, layout.pos.row=vpr))
    } else {
      pushViewport(dataViewport(xData=xlimdata, yData=ylim, extension=0,clip="off",
                                layout.pos.col=1, layout.pos.row=vpr))
      grid.yaxis(gp=gpar(cex=cexAxis))
      pushViewport(dataViewport(xData=xlimdata, yData=ylim, extension=0, clip="on",
                                layout.pos.col=1, layout.pos.row=vpr))
    }
  }
  
  defaultColors <- c("plus" = "#081d58", "minus" = "#081d58",
                     "duplicated" = "grey", "cp" = "#555555",
                     "ci" = "#777777", "highlight" = "red",
                     "threshold" = "grey")
  if (is.null(sampleColor))
    sampleColor <- defaultColors[strand]# "#081d58"
  
  ord  <- sort(which(dat$flag==0))# , which(dat$flag==0))
  colFlag <- ifelse(dat$flag==0, sampleColor, defaultColors["duplicated"])
  
### draw zero line
  if (zeroLine)
    grid.lines(y=unit(0, "native"), gp=gpar(col="#000000", lty=2, lwd=2))

  if (strand=="minus")
    dat$y <- -1*dat$y

  #grid.points(dat$x, y=dat$y,  pch=ipch, size=pointSize,gp=gpar(col=colFlag))
  if (length(dat$y)>0){ # do we have mapped reads on this strand?
  for (i in 1:length(dat$x.start))
    with(dat,
         grid.polygon(x=c(x.start[i],x.end[i], x.end[i],x.start[i]),
                      y=c(0, 0, y[i], y[i]), default.units="native",
                      gp=gpar(fill=colFlag[i])))#, col=colFlag[i] )))
    # grid.segments(x0=dat$x.start, x1=dat$x.end, y0=dat$y, y1=dat$y,
    #               default.units="native", gp=gpar(col=colFlag))
  }
  if(endVP)
     popViewport(2)
} ## plotReads
