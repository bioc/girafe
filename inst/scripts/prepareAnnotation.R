### Joern Toedling                                   September 2009

### This script shows how to create an object that holds information about
###  annotated genome features (positions of genes, ncRNAs etc.). Such an
###  object is required for some functionalities in package 'girafe'.

### In this case the constructed object will be of class
### 'Genome_intervals_stranded'

### annotation data from three sources is considered:
### 1. miRBase for microRNA sequences;
##     downloaded genomic coordinates of miRNA-precursors
###    URL  ftp://mirbase.org/pub/mirbase/13.0/genomes/mmu.gff
### 2. the RFAM database of annotated non-coding RNA families
###    URL  http://rfam.sanger.ac.uk
### 3. the MGI database of all mouse-related gene annotation
###    URL  http://www.informatics.jax.org/

### This script may be distributed according to the terms of
###  the Artistic License 2.0

### load the library:
library("girafe")


###-----------------------------------------------------------------
### I. ANNOTATION FROM miRBase
###-----------------------------------------------------------------
mm.mir <- readGff3("mmu_miR13.gff", isRightOpen=FALSE)
mm.mir$source <- rep("miRBase13", nrow(mm.mir))
mm.mir$family <- getGffAttribute(mm.mir, "ID")
mm.mir$class <- rep("miRNA", nrow(mm.mir))
mm.mir$"seq_name" <- factor(paste("chr", mm.mir$"seq_name",sep=""))

### add alias entry:
mm.mir$gffAttributes <- paste(mm.mir$gffAttributes, " Alias=", getGffAttribute(mm.mir, "ID"), sep="")


###-----------------------------------------------------------------
### II. ANNOTATION FROM RFAM DATABASE
###-----------------------------------------------------------------

### First we downloaded the annotation of ncRNAs in the mouse genome
###  from the RFAM database in form of the file
### ftp://ftp.sanger.ac.uk/pub/databases/Rfam/CURRENT/genome.gff3.tar.gz

### Files fro M_musculus were combined (via 'cat') into one file
###  containing the complete RFAM annotation for mouse

### Second we also downloaded a file that contains the description
###  of the RFAM families at
### ftp://ftp.sanger.ac.uk/pub/databases/Rfam/CURRENT/database_files/rfam.txt.gz

### read in the annotation file:
mm.rfam <- readGff3("RFAM_Mmusculus.gff", isRightOpen=FALSE)
### select only ncRNAs here:
mm.rfam <- mm.rfam[annotation(mm.rfam)$type=="ncRNA", drop=TRUE]

### read the family description file
rfam.desc <- read.delim("rfam.txt", header=FALSE)

### use this to replace the family identifiers in the annotation data.frame
###  by the actual names of the families:
mm.rfam$family <- getGffAttribute(mm.rfam, "Name")

rfam.names <- rfam.desc$V18
## in some cases there are problems in the file that the name is
##  in column 19 and not 18; column 18 only contains numbers in that case
problem.idx <- which(!is.na(as.numeric(rfam.names)))
if (length(problem.idx)>0)
  rfam.names[problem.idx] <-  rfam.desc$V19[problem.idx]
rfam.names <- gsub(";$","", rfam.names)
mm.rfam$class <- rfam.names[match(mm.rfam$family, rfam.desc$V3)]
### check whether all entries look reasonable
table(mm.rfam$class)

mm.rfam$class <- gsub("^Gene;","", mm.rfam$class)

### summarize the different classes of snoRNAs into one class
mm.rfam$class[grep("snoRNA", mm.rfam$class)] <- "snoRNA"

### now finally replace the missing family names (Y_RNA, 7SK RNA)
i <- which(mm.rfam$class=="Gene")
mm.rfam$class[i] <- getGffAttribute(mm.rfam, "Alias")[i]

#### remove miRNA entries (we have them from miRBase)
keep <- !mm.rfam$class %in% c("miRNA")
mm.rfam <- mm.rfam[keep]

## get an overview of which kinds of classes we have in the
###  annotation data and how often
table(mm.rfam$class)

### use the more common chromosome identifiers of the format chrX
###  instead of the RFAM-specific M_musculus_X
mm.rfam$"seq_name" <- factor(gsub("^M_musculus_","chr", mm.rfam$"seq_name"))
mm.rfam$class <- gsub(";splicing$", "", mm.rfam$class)


###-----------------------------------------------------------------
### III. GENE ANNOTATION FROM THE MGI DATABASE
###-----------------------------------------------------------------

### Second from the MGI database, we downloaded a GFF file containing
###  the genomic coordinates of all genes, pseudogenes in the mouse genome
### URL:  ftp://ftp.informatics.jax.org/pub/reports/MGI_GTGUP.gff
mm.mgi <- readGff3("MGI_GTGUP.gff", isRightOpen=FALSE)

### discard microRNA entries (we have them from miRBase) and gene models
keep <- !(annotation(mm.mgi)$type %in% c("microRNA","GeneModel"))
mm.mgi <- mm.mgi[keep]
### make sure both have the elements in their annotation
mm.mgi$family <- gsub("^.*;","", annotation(mm.mgi)$gffAttributes)
mm.mgi$class <- mm.mgi$type

### set standard value for phase '.' instead of NA
mm.mgi$phase[which(is.na(mm.mgi$phase))] <- "."

### some miRNAs had the unacceptable strand entry "None", but since
#    these are discarded there should not be a problem any longer
strand(mm.mgi) <- factor(as.character(strand(mm.mgi)), levels=c("-","+"))

stopifnot(validObject(mm.mgi))


###--------------------------------------------------------------
### IV. COMBINE THE OBJECTS INTO ONE ANNOTATION OBJECT
###--------------------------------------------------------------
mm.gi <- c(mm.mir, mm.rfam, mm.mgi)

### remove non-existant types
mm.gi$type <- factor(as.character(mm.gi$type))

### set score and phase
mm.gi$score <- rep(".", nrow(mm.gi))
mm.gi$phase <- rep(".", nrow(mm.gi))

### add entry Alias=
noAlias <- grep("Alias", mm.gi$gffAttributes, invert=TRUE)
mm.gi$gffAttributes[noAlias] <- gsub(";", ";Alias=", mm.gi$gffAttributes[noAlias])

### only use a restricted set of feature types for plotting:
allTypes <- as.character(mm.gi$class)
j <- which(allTypes %in% c("7SK", "GRIK4_3p_UTR", "NRON", "Telomerase-vert", "Vault", "Y_RNA"))
allTypes[j] <- "ncRNA"
j <- grep("gene", allTypes, ignore.case=TRUE)
allTypes[j] <- tolower(allTypes[j])
mm.gi$type <- allTypes

### get overview of types:
table(mm.gi$type)

### final check:
stopifnot(validObject(mm.gi))

### save the result object 'mm.gi':
save(mm.gi, file="anno_mm_genint.RData")



#####################################################################
### ALTERNATIVE USING PACKAGE rtracklayer    ########################
#####################################################################
# queries UCSC annotation tracks:

# 1. get positions of miRNA-precursors (as form miRBase)

session <- browserSession() # connect to UCSC
ucscGenomes() ## check which genome assemblies are available
genome(session) <- "mm9"
trackNames(session) ## check which annotation tracks are there
query <- ucscTableQuery(session, "miRNA")
# for seeing the result as a data.frame:
head(getTable(query))

## these data.frames could be converted to genome_intervals_stranded
##  objects and then used for annotation
