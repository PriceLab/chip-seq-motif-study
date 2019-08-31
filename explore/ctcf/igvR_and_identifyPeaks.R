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
if(!exists("findBindingSites"))
  source("findBindingSites.R")
#---------------------------
bamFile="sliceGSM749704_hg19_wgEncodeUwTfbsGm12878CtcfStdAlnRep1.bam"
tf="CTCF"
chrom= "chr19"
start=8993839
end=52184575
motifName="MA0139.1"
interestingSpot1<-"chr19:42,257,038-42,258,277"
interestingSpot2<-"chr19:43,033,314-43,034,146"
interestingSpot3<-"chr19:19,256,647-19,257,529"
interestingSpot4<-"chr19:42,210,396-42,211,437"
interestingSpot5<-"chr19:35,632,204-35,634,113"
interestingSpot6<-"chr19:42,937,212-42,945,062"

igv <- igvR() #connects to browser
setGenome(igv, "hg19")
showGenomicRegion(igv, interestingSpot1) 
idx <- indexBam(bamFile)

tbls<-identifyPeaks(bamFile, chrom, start, end, motifName) #gets my narrowpeak & broadpeak files 
tbl.narrow <- tbls[[1]]
tbl.broad <- tbls[[2]]
tbl.hits <- tbls[[3]]

tbl.narrow[,1]<-as.character(tbl.narrow[,1])
tbl.broad[,1]<-as.character(tbl.broad[,1])
narrowTrack <- DataFrameQuantitativeTrack("narrow", tbl.narrow[, c("chrom", "start.loc", "end.loc", "score")], color="red",autoscale=TRUE)
broadTrack <- DataFrameQuantitativeTrack("broad", tbl.broad[, c("chrom", "start.loc", "end.loc", "score")], color="green",autoscale=TRUE)
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
pileupTrack <- GenomicAlignmentTrack("ChIP-Seq PileUp", x)
displayTrack(igv, pileupTrack)

motifTrack <- DataFrameQuantitativeTrack("MA0139.1", tbl.hits[, c("chrom", "motifStart", "motifEnd", "motifRelativeScore")], color="blue",autoscale=TRUE)
displayTrack(igv, motifTrack)
