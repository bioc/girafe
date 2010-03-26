###################################################
### chunk number 1: prepare
###################################################
options(length=60, stringsAsFactors=FALSE)
set.seed(123)
options(SweaveHooks=list(
   along=function() par(mar=c(2.5,4.2,4,1.5), font.lab=2),
   pie=function() par(mar=c(0, 0, 0, 3.7), font=2)))


###################################################
### chunk number 2: loadpackage
###################################################
library("girafe")
library("RColorBrewer")


###################################################
### chunk number 3: setUp
###################################################
exDir <- system.file("extdata", package="girafe")
load(file.path(exDir, "anno_mm_genint.RData"))


###################################################
### chunk number 4: loadReads
###################################################
ra23.wa  <- readFastq(dirPath=exDir, pattern=
                      "aravinSRNA_23_plus_adapter_excerpt.fastq")


###################################################
### chunk number 5: showReads
###################################################
show(ra23.wa)


###################################################
### chunk number 6: trimAdapter
###################################################
adapter <- "CTGTAGGCACCATCAAT"
ra23.na  <- trimAdapter(ra23.wa, adapter)
show(ra23.na)


###################################################
### chunk number 7: readAligned
###################################################
exA   <- readAligned(dirPath=exDir, type="Bowtie", 
   pattern="aravinSRNA_23_no_adapter_excerpt_mm9_unmasked.bwtmap")
show(exA)


###################################################
### chunk number 8: convertAligned
###################################################
exAI <- as(exA, "AlignedGenomeIntervals")
organism(exAI) <- "Mm"


###################################################
### chunk number 9: showExAI
###################################################
show(exAI)


###################################################
### chunk number 10: tabChromosomes
###################################################
table(seq_name(exAI))


###################################################
### chunk number 11: showSubset
###################################################
detail(exAI[seq_name(exAI)=="chrMT"])


###################################################
### chunk number 12: setupToy
###################################################
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
detail(reduce(D))


###################################################
### chunk number 14: showReduceData
###################################################
S <- exAI[seq_name(exAI)=="chrX" & exAI@matches==1 & exAI[,1]>1e8]
detail(S)


###################################################
### chunk number 15: showReduceData2
###################################################
detail(reduce(S))


###################################################
### chunk number 16: plotAI
###################################################
plot(exAI, mm.gi, chr="chrX", start=50400000, end=50410000)


###################################################
### chunk number 17: examplePerWindow
###################################################
exPX  <- perWindow(exAI, chr="chrX", winsize=1e5, step=0.5e5)
head(exPX[order(exPX$n.overlap, decreasing=TRUE),])


###################################################
### chunk number 18: exportBed
###################################################
export(exAI, con=file.path(tempdir(), "export.bed"),
       format="bed", name="example_reads",
       description="Example reads",
       color="100,100,255", visibility="pack")


###################################################
### chunk number 19: showExport
###################################################
readLines(file.path(tempdir(), "export.bed"), n=4)


###################################################
### chunk number 20: getIntervalOverlap
###################################################
exOv <- interval_overlap(exAI, mm.gi)


###################################################
### chunk number 21: tableOverlap
###################################################
table(listLen(exOv))


###################################################
### chunk number 22: show12Elements
###################################################
getGffAttribute(mm.gi[exOv[[which(listLen(exOv)==12)]]],"Alias")[,1]


###################################################
### chunk number 23: computeTabOv
###################################################
(tabOv <- table(as.character(mm.gi$type)[unlist(exOv)]))


###################################################
### chunk number 24: displayPie
###################################################
my.cols <- brewer.pal(length(tabOv), "Set3")
pie(tabOv, col=my.cols, radius=0.95)


###################################################
### chunk number 25: sessionInfo
###################################################
toLatex(sessionInfo())


