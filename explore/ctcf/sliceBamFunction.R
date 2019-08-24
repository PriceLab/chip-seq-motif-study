library(Rsamtools)
library(RUnit)
library(org.Hs.eg.db)

#----------------------------------------------------------------------------------------------------
runTests <- function()
{
  test_sliceBamFile()
  
} # runTests
#----------------------------------------------------------------------------------------------------
sliceBamFile <- function(bamFile, chrom, start=NA, end=NA, newBamFileName)
{
  if(is.na(start)){
    start<- 1
  }
  
  if(is.na(end)){
    truncated.chromosome.name <- as.character(sub("chr", "", chrom))
    end <- org.Hs.egCHRLENGTHS[[truncated.chromosome.name]]
  }
  which <- GRanges(seqnames = chrom, ranges = IRanges(start, end))
  what <- c("rname", "strand", "pos", "qwidth", "seq")
  param <- ScanBamParam(which = which, what=what)
  
  bam <- scanBam(bamFile, param=param)
  
  filterBam(bamFile, newBamFileName, param=param)
  
} # sliceBamFile
#----------------------------------------------------------------------------------------------------
test_sliceBamFile<- function()
{
  message(sprintf("--- test_sliceBamFile"))
  
  bamFile <- "GSM749704_hg19_wgEncodeUwTfbsGm12878CtcfStdAlnRep1.bam"
  indexFile <- sprintf("%s.bai", bamFile)
  checkTrue(file.exists(bamFile))
  checkTrue(file.exists(indexFile))
  
  sliced.file <- "chr19-ctcf.bam"
  sliceBamFile(bamFile,"chr19", 1, 59128983, sliced.file)
  checkTrue(file.exists(sliced.file))
  
} # test_sliceBamFile
#----------------------------------------------------------------------------------------------------
if(!interactive())
  runTests()

