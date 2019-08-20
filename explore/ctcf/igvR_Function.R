runTests<- function()
{
  test_igvBam()
}

igvBam <- function(FileBam, tf, chromosome, narrowPeakFilefromMACS2, broadPeakFilefromMACS2)
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
  tbl.Np <- tables[[1]]
  tbl.Bp <- tables[[2]]
  colnames(tbl.Np) <- c("chrom", "start", "end", "name", "score","width")
  colnames(tbl.Bp) <- c("chrom", "start", "end", "name", "score", "width")
  
  
  tbl.Bp.Chrom <- tbl.bp[(grep(chromosome, tbl.Bp$chrom)[1]):((grep(chromosome, tbl.Bp$chrom)[1])+(length(grep(chromosome, tbl.Bp$chrom)))),]
  tbl.BestBroad.Chrom <-subset(tbl.Bp.Chrom, score == max(tbl.bp.chrom$score))
  
  
  start.loc.Specific= tbl.bestBroad.Chrom [1,'start']
  start.loc=  start.loc.Specific-5000
  end.loc.Specific= tbl.bestBroad.Chrom [1,'end']
  end.loc=end.loc.Specific+5000
    
  tbl.Bp.Chrom.Sub <- subset(tbl.bp.chrom, start.loc < start & end> end.loc)
  
  tbl.Np.Chrom <- tbl.Np[(grep(chromosome, tbl.Np$chrom)[1]): ((grep(chromosome, tbl.Np$chrom)[1])+(length(grep(chromosome, tbl.Np$chrom)))),]
  tbl.Np.Chrom.Sub <- subset(tbl.Np.Chrom,  start.loc < start & end> end.loc)

  NarrowTrack <- DataFrameAnnotationTrack("NarrowTrack", tbl.np.chrom.sub, color="red", displayMode="SQUISHED")
  BroadTrack <- DataFrameAnnotationTrack("BroadTrack", tbl.bp.chrom.sub, color="green", displayMode="SQUISHED")
  displayTrack(igv, NarrowTrack)
  displayTrack(igv, BroadTrack)
  
  
  
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
  
  igvBam(FileBam, tf, chromosome, narrowPeakFilefromMACS2, broadPeakFilefromMACS2)
  
  checkTrue(file.exists(FileBam))
  checkTrue(c("igvAnnotationTrack") %in% is(NarrowpTrack))
  checkTrue(c("igvAnnotationTrack") %in% is(BroadpTrack))
}
  
if(!interative()) #you can run tests from terminal
  runTests()




