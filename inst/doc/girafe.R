###################################################
### chunk number 1: prepare
###################################################
#line 65 "girafe.Rnw"
options(length=60, stringsAsFactors=FALSE)
set.seed(123)
options(SweaveHooks=list(
   along=function() par(mar=c(2.5,4.2,4,1.5), font.lab=2),
   pie=function() par(mar=c(0, 0, 0, 3.7), font=2)))


###################################################
### chunk number 2: loadpackage
###################################################
#line 88 "girafe.Rnw"
library("girafe")
library("RColorBrewer")


###################################################
### chunk number 3: setUp
###################################################
#line 102 "girafe.Rnw"
exDir <- system.file("extdata", package="girafe")
load(file.path(exDir, "anno_mm_genint.RData"))


###################################################
### chunk number 4: loadReads
###################################################
#line 111 "girafe.Rnw"
ra23.wa  <- readFastq(dirPath=exDir, pattern=
                      "aravinSRNA_23_plus_adapter_excerpt.fastq")


###################################################
### chunk number 5: showReads
###################################################
#line 116 "girafe.Rnw"
show(ra23.wa)


###################################################
### chunk number 6: trimAdapter
###################################################
#line 127 "girafe.Rnw"
adapter <- "CTGTAGGCACCATCAAT"
ra23.na  <- trimAdapter(ra23.wa, adapter)
show(ra23.na)


###################################################
### chunk number 7: readAligned
###################################################
#line 145 "girafe.Rnw"
exA   <- readAligned(dirPath=exDir, type="Bowtie", 
   pattern="aravinSRNA_23_no_adapter_excerpt_mm9_unmasked.bwtmap")
show(exA)


###################################################
### chunk number 8: convertAligned
###################################################
#line 155 "girafe.Rnw"
exAI <- as(exA, "AlignedGenomeIntervals")
organism(exAI) <- "Mm"


###################################################
### chunk number 9: showExAI
###################################################
#line 162 "girafe.Rnw"
show(exAI)


###################################################
### chunk number 10: tabChromosomes
###################################################
#line 167 "girafe.Rnw"
table(seq_name(exAI))


###################################################
### chunk number 11: showSubset
###################################################
#line 173 "girafe.Rnw"
detail(exAI[seq_name(exAI)=="chrMT"])


###################################################
### chunk number 12: setupToy
###################################################
#line 217 "girafe.Rnw"
D <- AlignedGenomeIntervals(
     start=c(1,3,4,5,8,10,11), end=c(5,5,6,8,9,11,13),
     chromosome=rep(c("chr1","chr2","chr3"), c(2,2,3)),
     strand=c("-","-","+","+","+","+","+"),
     sequence=c("ACATT","ACA","CGT","GTAA","AG","CT","TTT"),
     reads=rep(1,7), matches=c(rep(1,6),3))

detail(D)


###################################################
### chunk number 13: showReduceToy
###################################################
#line 229 "girafe.Rnw"
detail(reduce(D))


###################################################
### chunk number 14: showReduceData
###################################################
#line 238 "girafe.Rnw"
S <- exAI[seq_name(exAI)=="chrX" & exAI@matches==1 & exAI[,1]>1e8]
detail(S)


###################################################
### chunk number 15: showReduceData2
###################################################
#line 246 "girafe.Rnw"
detail(reduce(S))


###################################################
### chunk number 16: reduceExample3
###################################################
#line 262 "girafe.Rnw"
S2 <- exAI[seq_name(exAI)=="chr11" & exAI@matches==1 & exAI[,1]>8e7]
detail(S2)
detail(reduce(S2, exact=TRUE))


###################################################
### chunk number 17: plotAI
###################################################
#line 281 "girafe.Rnw"
plot(exAI, mm.gi, chr="chrX", start=50400000, 
     end=50410000, show="minus")


###################################################
### chunk number 18: examplePerWindow
###################################################
#line 320 "girafe.Rnw"
exPX  <- perWindow(exAI, chr="chrX", winsize=1e5, step=0.5e5)
head(exPX[order(exPX$n.overlap, decreasing=TRUE),])


###################################################
### chunk number 19: exportBed eval=FALSE
###################################################
## #line 341 "girafe.Rnw"
## export(exAI, con="export.bed",
##        format="bed", name="example_reads",
##        description="Example reads",
##        color="100,100,255", visibility="pack")


###################################################
### chunk number 20: getIntervalOverlap
###################################################
#line 374 "girafe.Rnw"
exOv <- interval_overlap(exAI, mm.gi)


###################################################
### chunk number 21: tableOverlap
###################################################
#line 380 "girafe.Rnw"
table(listLen(exOv))


###################################################
### chunk number 22: show12Elements
###################################################
#line 386 "girafe.Rnw"
getGffAttribute(mm.gi[exOv[[which(listLen(exOv)==12)]]],"Alias")[,1]


###################################################
### chunk number 23: computeTabOv
###################################################
#line 393 "girafe.Rnw"
(tabOv <- table(as.character(mm.gi$type)[unlist(exOv)]))


###################################################
### chunk number 24: displayPie
###################################################
#line 399 "girafe.Rnw"
my.cols <- brewer.pal(length(tabOv), "Set3")
pie(tabOv, col=my.cols, radius=0.95)


###################################################
### chunk number 25: multicoreShow eval=FALSE
###################################################
## #line 508 "girafe.Rnw"
## library("multicore")
## options("cores"=4) # adjust to your machine
## covAI <- coverage(exAI, byStrand=TRUE)


###################################################
### chunk number 26: sessionInfo
###################################################
#line 556 "girafe.Rnw"
toLatex(sessionInfo(), locale=FALSE)


