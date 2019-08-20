runTests<- function()
{
  test_igvBam()
}

igvBam <- function(FileBam, tf, chromsome, narrowPeakFilefromMACS2, broadPeakFilefromMACS2)
{

  library(RUnit)
  library(igvR)
  library(VariantAnnotation)
  
  
  idx = indexBam(FileBam)
  igv <- igvR() #connects to browser
  setGenome(igv, "hg19")
  showGenomicRegion(igv, tf)
  
  unzippedTables <- lapply (c(narrowPeakFilefromMACS2,broadPeakFilefromMACS2), read.table)
  tables <- lapply(unzippedTables,as.data.frame)
  tbl.np <- tables[[1]]
  tbl.bp <- tables[[2]]
  colnames(tbl.np) <- c("chrom", "start", "end", "name", "score","width")
  colnames(tbl.bp) <- c("chrom", "start", "end", "name", "score", "width")
  
  
  tbl.bp.chrom <- tbl.bp[(grep(chromosome, tbl.bp$chrom)[1]):((grep(chromosome, tbl.bp$chrom)[1])+(length(grep(chromosome, tbl.bp$chrom)))),]
  tbl.bestBroad.Chrom <-subset(tbl.bp.chrom, score == max(tbl.bp.chrom$score))
  
  
  start.loc.Specific= tbl.bestBroad.Chrom [1,'start']
  start.loc=  start.loc.Specific-5000
  end.loc.Specific= tbl.bestBroad.Chrom [1,'end']
  end.loc=end.loc.Specific+5000
    
  tbl.bp.chrom.sub <- subset(tbl.bp.chrom, start.loc < start & end> end.loc)
  
  tbl.np.chrom <- tbl.np[(grep(chromosome, tbl.np$chrom)[1]): ((grep(chromosome, tbl.np$chrom)[1])+(length(grep(chromosome, tbl.np$chrom)))),]
  tbl.np.chrom.sub <- subset(tbl.np.chrom,  start.loc < start & end> end.loc)

  NarrowpTrack <- DataFrameAnnotationTrack("NarrowpTrack", tbl.np.chrom.sub, color="red", displayMode="SQUISHED")
  BroadpTrack <- DataFrameAnnotationTrack("BroadpTrack", tbl.bp.chrom.sub, color="green", displayMode="SQUISHED")
  displayTrack(igv, NarrowpTrack)
  displayTrack(igv, BroadpTrack)
  
  
  
  which <- GRanges(seqnames = chromosome, ranges = IRanges(start.loc, end.loc))
  param <- ScanBamParam(which=which)
  
  x <- readGAlignments(FileBam, use.names=TRUE, param=param)
  track <- GenomicAlignmentTrack("DNAse", x)
  displayTrack(igv, track)
}

#test
test_igvBam <- function()
{
  message(sprintf("--- test_igvBam"))
  
  FileBam="GSM749704_hg19_wgEncodeUwTfbsGm12878CtcfStdAlnRep1.bam"
  tf="CTCF"
  chromosome= "chr11"
  narrowPeakFilefromMACS2= "ctcf__peaks.narrowPeak"
  broadPeakFilefromMACS2= "ctcf__peaks.broadPeak"
  
  igvBam(FileBam, tf, chromsome, narrowPeakFilefromMACS2, broadPeakFilefromMACS2)
  
  checkTrue(file.exists(FileBam))
  checkTrue(c("igvAnnotationTrack") %in% is(NarrowpTrack))
  checkTrue(c("igvAnnotationTrack") %in% is(BroadpTrack))
}
  
if(!interative())
  runTests()




