library(RPostgreSQL)
library(FimoClient)

pfms <- as.list(query(MotifDb, c("sapiens", "TEAD1", "jaspar2018", "MA0090.2")))
length(pfms)
fc <- FimoClient("khaleesi", 60010)
db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="hg38", host="khaleesi")

tf <- "TEAD1"
chrom <- "chr19"

#-----------------------------------------------------------------------
# just 3 TEAD1 ChIP-seq hits in this 14.4k region, good for a quick test
#-----------------------------------------------------------------------
start.loc <- 44902547
end.loc   <- 44917011

#-----------------------------------------------------------------------
# 2030 TEAD1 ChIP-seq hits in this 58.6M region on chr19
#-----------------------------------------------------------------------
start.loc <- 1
end.loc <- 58617616

query <- sprintf("select * from chipseq where tf='%s' AND chrom='%s' AND start >= %d AND endpos <= %d", tf, chrom, start.loc, end.loc)
tbl.tf <- dbGetQuery(db, query)
dim(tbl.tf)

tbl.tf.adjusted <- tbl.tf[, c("chrom", "peakStart", "peakEnd")]
colnames(tbl.tf.adjusted) <- c("chrom", "start", "end")
tbl.tf.adjusted$start <- tbl.tf.adjusted$start - 9
tbl.tf.adjusted$end <- tbl.tf.adjusted$end + 9
dim(tbl.tf.adjusted)
tbl.motifs <- requestMatchForRegions(fc, tbl.tf.adjusted, "hg38", .05)
dim(tbl.motifs)

#------------------------------------------------------------
# now associate these fimo hits with the chip hits
#------------------------------------------------------------
gr.chip <- GRanges(tbl.tf.adjusted)
gr.motifs <- GRanges(tbl.motifs)
tbl.overlaps <- as.data.frame(findOverlaps(gr.motifs, gr.chip, type="any"))
colnames(tbl.overlaps) <- c("motif", "chip")
chips.matched <- length(unique(tbl.overlaps$chip))
chips.unmatched <- nrow(tbl.tf.adjusted) - chips.matched
motifs.matched <- length(unique(tbl.overlaps$motif))
fimo.pValScores <- -log10(tbl.motifs$pValue)
dist <- fivenum(fimo.pValScores)
printf("found motif match for %d/%d ChIP hits (%5.2f%%)",
       chips.matched, nrow(tbl.tf.adjusted), 100 * chips.matched/nrow(tbl.tf.adjusted))
printf("-log10(fimo pValue) range of %f %f %f %f %f",
       dist[1], dist[2], dist[3], dist[4], dist[5])
fimo.pValScores.withFailures <- c(fimo.pValScores, rep(-1, chips.unmatched))
hist(fimo.pValScores.withFailures,
     main=sprintf("JASPAR TEAD1 FIMO motif match to ChIP hits\n with %d failures (chr19)", chips.unmatched),
     xlab="-log10(FIMO pValue)")

