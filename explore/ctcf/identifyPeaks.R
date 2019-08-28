# identifyPeaks.R:  call dockerized macs2, on the specified slice of a bamfile, return a list of narrow and broad peaks
#------------------------------------------------------------------------------------------------------------------------
library(Rsamtools)
library(RUnit)
library(org.Hs.eg.db)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("sliceBamFile"))
  source("sliceBamFile.R")
#------------------------------------------------------------------------------------------------------------------------
identifyPeaks<- function(bamFile, chrom, start, end)
{
  
  truncated.bamFile.name <- as.character(sub(".bam","", bamFile))
  
  broadpeaks.cmd <- sprintf("docker run -v $PWD:/data/ fooliu/macs2 callpeak -t /data/%s -c --broad -g hs --broad-cutoff 0.1 -n %s --outdir /data/",
                            bamFile,truncated.bamFile.name)
  
  system(broadpeaks.cmd)
  narrowpeaks.cmd <- sprintf("docker run -v $PWD:/data/ fooliu/macs2 callpeak -t /data/%s -c  -n %s --outdir /data/",
                             bamFile,truncated.bamFile.name)
  system(narrowpeaks.cmd)
  
  unzippedTables<- lapply(sprintf(c("%s_peaks.narrowPeak","%s_peaks.broadPeak"),truncated.bamFile.name), read.table)
  tables <- lapply(unzippedTables,as.data.frame)
  tbl.narrow <- tables[[1]]
  tbl.broad <- tables[[2]]
  colnames(tbl.narrow) <- c("chrom", "start.loc", "end.loc", "name", "score","width")
  colnames(tbl.broad) <- c("chrom", "start.loc", "end.loc", "name", "score", "width")
  
  tbls<- list(tbl.narrow, tbl.broad)
  names(tbls)<- c("narrow","broad")
  return(tbls)
  
} # identifyPeaks
#------------------

# test_identifyPeaks.R
source("identifyPeaks.R")
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
  test_identifyPeaks()
  
} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_identifyPeaks <- function()
{
  bamFile <- "sliceGSM749704_hg19_wgEncodeUwTfbsGm12878CtcfStdAlnRep1.bam"
  
  tbls <- identifyPeaks(bamFile, chrom="chr19", start=53160025, end=56170311)
  checkTrue(is.list(tbls))
  checkTrue(all(c("narrow", "broad") %in% names(tbls)))
  
  tbl.narrow <- tbls$narrow
  tbl.broad <- tbls$broad
  
  checkTrue(is.data.frame(tbl.narrow))
  checkTrue(is.data.frame(tbl.broad))
  
  checkTrue(nrow(tbl.narrow) > 0)
  checkTrue(nrow(tbl.broad) > 0)
  
  # many more checks to come
  
} # test_identifyPeaks
#-----------------------
