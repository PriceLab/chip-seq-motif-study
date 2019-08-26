identifyPeaks<- function(bamFile, chrom, start.loc, end.loc, newBamFileName)
{
  if (exists("sliceBamFile")){
    sliceBamFile(bamFile, chrom, start.loc, end.loc, newBamFileName)
    } else {
      source("sliceBamFile.R")
      sliceBamFile(bamFile, chrom, start.loc, end.loc, newBamFileName)
      }
  truncated.new.bamFile.name <- as.character(sub(".bam","", newBamFileName))
  system(sprintf("docker run -v $PWD:/data/ fooliu/macs2 callpeak -t /data/%s -c --broad -g hs --broad-cutoff 0.1 -n %s --outdir /data/", newBamFileName,truncated.new.bamFile.name))
  system(sprintf("docker run -v $PWD:/data/ fooliu/macs2 callpeak -t /data/%s -c  -n %s --outdir /data/", newBamFileName,truncated.new.bamFile.name))
  
  unzippedTables<- lapply(sprintf(c("%s_peaks.narrowPeak","%s_peaks.broadPeak"),truncated.new.bamFile.name), read.table)
  tables <- lapply(unzippedTables,as.data.frame)
  tbl.narrow <- tables[[1]]
  tbl.broad <- tables[[2]]
  colnames(tbl.narrow) <- c("chrom", "start", "end", "name", "score","width")
  colnames(tbl.broad) <- c("chrom", "start", "end", "name", "score", "width")
  tbl.broad.chrom <- tbl.broad[(grep(chrom, tbl.broad$chrom)[1]):((grep(chrom, tbl.broad$chrom)[1])+(length(grep(chrom, tbl.broad$chrom)))),]
  tbl.broad.chrom.sub <- subset(tbl.broad.chrom, start.loc < start & end> end.loc)
  
  tbl.narrow.chrom <- tbl.narrow[(grep(chrom, tbl.narrow$chrom)[1]): ((grep(chrom, tbl.narrow$chrom)[1])+(length(grep(chrom, tbl.narrow$chrom)))),]
  tbl.narrow.chrom.sub <- subset(tbl.narrow.chrom,  start.loc < start & end> end.loc)
  
  if(is.na(end.loc)){
    tbls<- list(tbl.narrow.chrom,tbl.broad.chrom)
    names(tbls)<- c("narrow","broad")
    return(tbls) 
    } else { 
      tbls<- list(tbl.narrow.chrom.sub,tbl.broad.chrom.sub)
      names(tbls)<- c("narrow","broad")
      return(tbls) 
    }
  }
