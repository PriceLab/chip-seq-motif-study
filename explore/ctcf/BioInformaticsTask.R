retrieveMatches <- function(tf, motifName, chromosome, start.loc=NA, end.loc=NA, displayHistogram=FALSE, peakMode="narrow", fimoThreshold)
{quietly=FALSE
quiet=FALSE
library(RPostgreSQL)
library(FimoClient)
library(MotifDb)
library(org.Hs.eg.db)
pfms <- as.list(query(MotifDb, c("sapiens", tf, "jaspar2018", motifName)))
length(pfms)
fc <- FimoClient("khaleesi", 60010, quiet=FALSE)
db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="hg38", host="khaleesi")
tf <- tf
chrom <- chromosome
#default to small and make 

if(is.na(start.loc)){
  start.loc=1
}

if(is.na(end.loc)){
  truncated.chromosome.name <- as.character(sub("chr", "", chromosome))
  end.loc <- org.Hs.egCHRLENGTHS[[truncated.chromosome.name]]
}

query <- sprintf("select * from chipseq where tf='%s' AND chrom='%s' AND start >= %d AND endpos <= %d", tf, chromosome, start.loc, end.loc)
tbl.tf <- dbGetQuery(db, query)
dim(tbl.tf)
colnames(tbl.tf) #what are the colnames 
print(colnames(tbl.tf))

if(peakMode=="narrow"){
  startColumn <- "peakStart"
  endColumn <- "peakEnd"
  tbl.tf["chrom", startColumn, endColumn]  
}
if(peakMode=="broad"){
  startColumn <- "start"
  endColumn <- "endpos"
  tbl.tf["chrom", startColumn, endColumn]  
}

tbl.tf.adjusted <- tbl.tf[, c("chrom", "peakStart", "peakEnd")]
colnames(tbl.tf.adjusted) <- c("chrom", "start", "end")
tbl.tf.adjusted$start <- tbl.tf[,"start"]
tbl.tf.adjusted$end <- tbl.tf[,"endpos"]
dim(tbl.tf.adjusted)

tbl.motifs <- requestMatchForRegions(fc, tbl.tf, "hg38", fimoThreshold)
dim(tbl.motifs)

gr.chip <- GRanges(tbl.tf.adjusted)
gr.motifs <- GRanges(tbl.motifs)
tbl.overlaps <- as.data.frame(findOverlaps(gr.motifs, gr.chip, type="any"))
colnames(tbl.overlaps) <- c("motif", "chip")
chips.matched <- length(unique(tbl.overlaps$chip))
chips.unmatched <- nrow(tbl.tf.adjusted) - chips.matched
motifs.matched <- length(unique(tbl.overlaps$motif))
fimo.pValScores <- -log10(tbl.motifs$pValue)
dist <- fivenum(fimo.pValScores)
printf <- function(...) print(noquote(sprintf(...)))
printf("found motif match for %d/%d ChIP hits (%5.2f%%)",
       chips.matched, nrow(tbl.tf.adjusted), 100 * chips.matched/nrow(tbl.tf.adjusted))
printf("-log10(fimo pValue) range of %f %f %f %f %f",
       dist[1], dist[2], dist[3], dist[4], dist[5])
fimo.pValScores.withFailures <- c(fimo.pValScores, rep(-1, chips.unmatched))
if (displayHistogram) {
  hist(fimo.pValScores.withFailures,
       main=sprintf("JASPAR %s FIMO motif match to ChIP hits\n with %d failures (chr19)",tf, chips.unmatched),
       xlab="-log10(FIMO pValue)") 
}
}
tbl.broadPeaks <- read.table(gzfile("GSM1010734_hg19_wgEncodeHaibTfbsMcf7CtcfcV0422111PkRep1.broadPeak"), header= TRUE) #reading the data we have
tbl.broadPeaks 
colnames(tbl.broadPeaks) <- c("chrom", "start", "end", "name", "score") #we are getting the columns that we want from the tbl
fivenum(tbl.broadPeaks$score) #summarizes the table in five numbers (minimum,lower-hinge, median, upper-hinge, maximum)
tbl.best <- subset(tbl.broadPeaks, score == max(tbl.broadPeaks$score)) #Find the highest scoring ChIP-seq hit reported in the data
tbl.best 

retrieveMatches("CTCF", "MA0139.1", "chr15", start.loc=39109648, end.loc=39110469, displayHistogram=FALSE, peakMode="broad", .05)
retrieveMatches("CTCF", "MA0139.1", "chr3", start.loc=184278744, end.loc=184279631, displayHistogram=FALSE, peakMode="broad", .05)
#Find a match of the CTCF motif to the DNA sequence in that region

#RESULTS: 
http://meme-suite.org/opal-jobs/appFIMO_5.0.41565285676197436389089/fimo.html
http://meme-suite.org/opal-jobs/appFIMO_5.0.41565285889843-1673161087/fimo.html


#Notes 
# - Marjorie Brand Laboratory- located in Canada, studies the regulation of gene expression at the level of transcription and chromatin
# - GEO (gene expression omnibus)- a database that has experimental data archived (primarily used as data storage)
# - MACS2 broadpeak format- a format of data 
# - hg19- a human reference genome (it's version 19)
# - GSM1010734- one of the ChIP-seq experiments archived on GEO 
# - file extension “.gz”   unzip tool “gunzip" - used to unzip .gz files, so that it can be read and used by R 

# hg19 sequence:  https://www.biostars.org/p/56/
# FIMO on the web: http://meme-suite.org/meme_5.0.4/tools/fimo
