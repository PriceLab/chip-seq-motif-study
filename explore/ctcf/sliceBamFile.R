# identifyPeaks.R:  call dockerized macs2, on the specified slice of a bamfile, return a list of narrow and broad peaks
#------------------------------------------------------------------------------------------------------------------------
library(Rsamtools)
library(RUnit)
library(org.Hs.eg.db)
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
