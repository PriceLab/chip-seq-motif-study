#simple script using identifyPeaks
#----------------------------------
library(RUnit)
library(igvR)
library(VariantAnnotation)
library(Rsamtools)
library(RUnit)
library(org.Hs.eg.db)
#---------------------------
if(!exists("identifyPeaks"))
  source("identifyPeaks.R")
#---------------------------
bamFile="sliceGSM749704_hg19_wgEncodeUwTfbsGm12878CtcfStdAlnRep1.bam"
tf="CTCF"
chrom= "chr19"
start=53160025
end=56170311

igv <- igvR() #connects to browser
setGenome(igv, "hg19")
showGenomicRegion(igv, tf) 
idx <- indexBam(bamFile)

tbls<-identifyPeaks(bamFile, chrom, start, end) #gets my narrowpeak & broadpeak files 
tbl.narrow <- tbls[[1]]
tbl.broad <- tbls[[2]]
narrowTrack <- DataFrameAnnotationTrack("narrowTrack", tbl.narrow, color="red", displayMode="SQUISHED")
broadTrack <- DataFrameAnnotationTrack("broadTrack", tbl.broad, color="green", displayMode="SQUISHED")
displayTrack(igv, narrowTrack)
displayTrack(igv, broadTrack)

start.loc=NA
end.loc=NA 
if(is.na(start.loc)){
  start.loc<- 1
  }

if(is.na(end.loc)){
  truncated.chromosome.name <- as.character(sub("chr", "", chrom))
  end.loc <- org.Hs.egCHRLENGTHS[[truncated.chromosome.name]]
  }

which <- GRanges(seqnames = chrom, ranges = IRanges(start.loc, end.loc))
param <- ScanBamParam(which=which)
x <- readGAlignments(bamFile, use.names=TRUE, param=param)
track <- GenomicAlignmentTrack("DNAse", x)
displayTrack(igv, track)
