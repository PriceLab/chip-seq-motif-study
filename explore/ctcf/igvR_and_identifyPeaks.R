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
start=53160025
end=56170311
motifName="MA0139.1"

igv <- igvR() #connects to browser
setGenome(igv, "hg19")
showGenomicRegion(igv, "chr19:1-59,128,983") 
idx <- indexBam(bamFile)

tbls<-identifyPeaks(bamFile, chrom, start, end, motifName) #gets my narrowpeak & broadpeak files 
tbl.narrow <- tbls[[1]]
tbl.broad <- tbls[[2]]
tbl.hits <- tbls[[3]]

narrowTrack <- DataFrameAnnotationTrack("narrow", tbl.narrow, color="red", displayMode="SQUISHED")
broadTrack <- DataFrameAnnotationTrack("broad", tbl.broad, color="green", displayMode="SQUISHED")
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

motifTrack <- DataFrameQuantitativeTrack("MA0139.1", tbl.hits[, c("chrom", "motifStart", "motifEnd", "motifRelativeScore")], color="blue",autoscale=TRUE,)
displayTrack(igv, motifTrack)
