sliceBamFile <- function(existingBamFileNam, chrom, start, end, newBamFileName)
{
  library(Rsamtools)
  which <- GRanges(seqnames = chrom, ranges = IRanges(start, end))
  what <- c("rname", "strand", "pos", "qwidth", "seq")
  param <- ScanBamParam(which = which, what=what)
  bamFile<- sprintf("~/github/igvR/inst/extdata/%s", existingBamFileNam)
  bam <- scanBam(bamFile, param=param)
  filterBam(bamFile, newBamFileName, param=param)
}

#unit test
testSliceBamFile<- function()
{
  message(sprintf("--- test_igvBam"))
  
  sliceBamFile ("GSM749704_hg19_wgEncodeUwTfbsGm12878CtcfStdAlnRep1.bam","chr19", 1, 59128983, "temp1.bam")
  checkTrue(file.exists("GSM749704_hg19_wgEncodeUwTfbsGm12878CtcfStdAlnRep1.bam"))
  checkTrue(file.exists("temp1.bam"))
  }
