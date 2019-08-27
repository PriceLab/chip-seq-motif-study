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

#  if (exists("sliceBamFile")){
#    sliceBamFile(bamFile, chrom, start, end, newBamFileName)
#  } else {
#    source("sliceBamFile.R")
#    sliceBamFile(bamFile, chrom, start, end, newBamFileName)
#  }

  browser()
  newBamFileName<-sprintf(("slice%s"), bamFile)

  sliceBamFile(bamFile, chrom, start, end, newBamFileName)
  truncated.new.bamFile.name <- as.character(sub(".bam","", newBamFileName))

  broadpeaks.cmd <- sprintf("docker run -v $PWD:/data/ fooliu/macs2 callpeak -t /data/%s -c --broad -g hs --broad-cutoff 0.1 -n %s --outdir /data/",
                            newBamFileName,truncated.new.bamFile.name)

  system(broadpeaks.cmd)
  narrowpeaks.cmd <- sprintf("docker run -v $PWD:/data/ fooliu/macs2 callpeak -t /data/%s -c  -n %s --outdir /data/",
                             newBamFileName,truncated.new.bamFile.name)
  system(narrowpeaks.cmd)

  unzippedTables<- lapply(sprintf(c("%s_peaks.narrowPeak","%s_peaks.broadPeak"),truncated.new.bamFile.name), read.table)
  tables <- lapply(unzippedTables,as.data.frame)
  tbl.narrow <- tables[[1]]
  tbl.broad <- tables[[2]]
  colnames(tbl.narrow) <- c("chrom", "start.loc", "end.loc", "name", "score","width")
  colnames(tbl.broad) <- c("chrom", "start.loc", "end.loc", "name", "score", "width")

  tbls<- list(tbl.narrow, tbl.broad)
  names(tbls)<- c("narrow","broad")
  return(tbls)

} # identifyPeaks
#------------------------------------------------------------------------------------------------------------------------
