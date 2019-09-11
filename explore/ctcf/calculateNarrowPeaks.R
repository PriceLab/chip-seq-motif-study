library(Rsamtools)
library(org.Hs.eg.db)
library(ChIPseqMotifMatch)
library(RUnit)
source("sliceBamFile.R")
#------------
test.bamFilename <- system.file(package="ChIPseqMotifMatch", "extdata", "GSM749704_hg19_Ctcf_chr19.bam")
idx <- indexBam(test.bamFilename)

slice.bamFile<- sliceBamFile(bamFile= test.bamFilename, chrom="chr19", start=4225703	, end=50258277, "sliced.bam")

truncated.slice.bamFile.name <- as.character(sub(".bam","", slice.bamFile))

narrowpeaks.cmd <- sprintf("docker run -v $PWD:/data/ fooliu/macs2 callpeak -t /data/%s -c  -n %s --outdir /data/",
                           slice.bamFile,truncated.slice.bamFile.name)
system(narrowpeaks.cmd)

unzippedTable<- read.table(gsub("%s",truncated.slice.bamFile.name,"%s_peaks.narrowPeak"))

tbl.narrow  <- as.data.frame(unzippedTable)

colnames(tbl.narrow) <- c("chrom", "start.loc", "end.loc", "name", "score","width", "foldChange", "pValue","qValue", "summitPosition")

head(tbl.narrow)
