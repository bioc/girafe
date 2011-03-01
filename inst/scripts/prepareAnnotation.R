### Joern Toedling                                   February 2011

### This script shows how to create an object that holds information about
###  annotated genome features (positions of genes, ncRNAs etc.). Such an
###  object is required for some functionalities in package 'girafe'.

### In this case the constructed object will be of class
### 'Genome_intervals_stranded'

### The table of annotated genome features is retrieved from
### the MGI database
###    URL  http://www.informatics.jax.org/

### This script may be distributed according to the terms of
###  the Artistic License 2.0

###-----------------------------------------------------------------
### 0. RETRIEVE ANNOTATION TABLE FROM MGI DATABASE
###-----------------------------------------------------------------
# load the library:
library("girafe")

###-----------------------------------------------------------------
### I. RETRIEVE ANNOTATION TABLE FROM MGI DATABASE
###-----------------------------------------------------------------

# Use "Genes and Markers Query Form" at MGI web site:
# http://www.informatics.jax.org/
# Select 
# 1. protein-coding genes, non-coding RNA genes and pseudogenes
# 2. all chromosomes
# 3. sort output by genome coordinates
# 4. Output: tab-delimited to FTP
# and hit "Search" to start assembly of the table.
# The resulting tab-delimited file can be then downloaded from
# the MGI FTP server.
# In this case, it is called "MGI_features_mm9_Feb2011.txt"

###-----------------------------------------------------------------
### II. IMPORT & POST-PROCESS TABLE IN R
###-----------------------------------------------------------------

mgi <- read.delim("MGI_features_mm9_Feb2011.txt")
table(mgi$"Feature.Type")

### remove unmapped features:
sum(is.na(mgi$start)) # how many to remove?
mgi <- mgi[!is.na(mgi$start),]

### consistency check: are all end sites >= start sites
###  (necessary, we noticed an error for MiR-1902 and had to correct
###   that in the file according to the miRBase annotation)
all(mgi$end >= mgi$start)

### rename MGI feature types 
mgi.types <- gsub(" gene$", "", mgi$"Feature.Type")
mgi.types <- gsub("protein coding", "gene", mgi.types)
table(mgi.types)

## which are the target types to consider:
accepted.types <- c("gene", "pseudogene", "miRNA", "lincRNA", "rRNA",
                    "snRNA", "snoRNA", "tRNA", "ncRNA")
# 'ncRNA' is used as a summary term for the
#  other small classes of RNA genes:
mgi.types <- factor(mgi.types, levels=accepted.types)
mgi.types[is.na(mgi.types)] <- "ncRNA"

## after processing we have the following feature types and numbers:
table(mgi.types)


###-----------------------------------------------------------------
### III. create AlignedGenomeIntervals object
###-----------------------------------------------------------------

mgi.gi <-
  GenomeIntervals(chromosome=paste('chr', mgi$Chr, sep=''),
                  start=mgi$start,
                  end=mgi$end,
                  strand=mgi$strand,
                  source=rep("MGI", nrow(mgi)),
                  type=as.character(mgi.types),
                  gffAttributes=paste('ACC="', mgi$MGI.ID, '"; ',
                    'Name="', mgi$Name, '"', sep=''),
                  ID=mgi$Symbol )

save(mgi.gi, file="mgi_gi.RData")
