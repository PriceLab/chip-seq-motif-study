runTests<- function()
{
  testIgvBam()
}

igvBam <- function(fileBam, tf, chromosome, narrowPeakFilefromMACS2, broadPeakFilefromMACS2)
{

  library(RUnit)
  library(igvR)
  library(VariantAnnotation)
  
  
  idx = indexBam(fileBam)
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
  tbl.bestBroad.chrom <-subset(tbl.bp.chrom, score == max(tbl.bp.chrom$score))
  
  
  start.loc.specific= tbl.bestBroad.chrom [1,'start']
  start.loc= start.loc.specific-5000
  end.loc.specific= tbl.bestBroad.chrom [1,'end']
  end.loc=end.loc.specific+5000
    
  tbl.bp.chrom.sub <- subset(tbl.bp.chrom, start.loc < start & end> end.loc)
  
  tbl.np.chrom <- tbl.np[(grep(chromosome, tbl.np$chrom)[1]): ((grep(chromosome, tbl.np$chrom)[1])+(length(grep(chromosome, tbl.np$chrom)))),]
  tbl.np.chrom.sub <- subset(tbl.np.chrom,  start.loc < start & end> end.loc)

  narrowTrack <- DataFrameAnnotationTrack("narrowTrack", tbl.np.chrom.sub, color="red", displayMode="SQUISHED")
  broadTrack <- DataFrameAnnotationTrack("broadTrack", tbl.bp.chrom.sub, color="green", displayMode="SQUISHED")
  displayTrack(igv, narrowTrack)
  displayTrack(igv, broadTrack)
  
  
  
  which <- GRanges(seqnames = chromosome, ranges = IRanges(start.loc, end.loc))
  param <- ScanBamParam(which=which)
  
  x <- readGAlignments(fileBam, use.names=TRUE, param=param)
  track <- GenomicAlignmentTrack("DNAse", x)
  displayTrack(igv, track)
}

#test
testIgvBam <- function()
{
  message(sprintf("--- test_igvBam"))
  
  fileBam="GSM749704_hg19_wgEncodeUwTfbsGm12878CtcfStdAlnRep1.bam"
  checkTrue(file.exists(fileBam))
  
  tf="CTCF"
  chromosome= "chr11"
  narrowPeakFilefromMACS2= "ctcf__peaks.narrowPeak"
  broadPeakFilefromMACS2= "ctcf__peaks.broadPeak"
  
  igvBam(fileBam, tf, chromosome, narrowPeakFilefromMACS2, broadPeakFilefromMACS2)
  
  checkTrue(c("igvAnnotationTrack") %in% is(narrowTrack))
  checkTrue(c("igvAnnotationTrack") %in% is(broadTrack))
}
  
if(!interactive()) #you can run tests from terminal
  runTests()




