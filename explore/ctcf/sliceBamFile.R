# sliceBamFile.R: extract a subset of a bam file
# developed originally to create a small yet representative bam file for the
# empirical motif-matching/ChIP-seq project
#----------------------------------------------------------------------------------------------------
library(Rsamtools)
library(RUnit)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_sliceBamFile()

} # runTests
#----------------------------------------------------------------------------------------------------
sliceBamFile <- function(bamFile, chrom, start, end, newBamFileName)
{
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
