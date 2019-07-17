library(RPostgreSQL)
library(FimoClient)
library(MotifDb)
library(org.Hs.eg.db)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_CTCF_2kbRegionOnChr19()

} # runTests
#----------------------------------------------------------------------------------------------------
retrieveMatches <- function(tf, motifName, chromosome, start.loc=NA, end.loc=NA, displayHistogram=FALSE,
                            fimoThreshold=1e-4, quiet=TRUE)
{
   pfms <- as.list(query(MotifDb, c("sapiens", tf, "jaspar2018", motifName)))
   length(pfms)
   fc <- FimoClient("khaleesi", 60010, quiet=FALSE)
   db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="hg38", host="khaleesi")

       # if not specified, assign chrom start and end - the entire lengt of the requested chromosome

   if(is.na(start.loc))
      start.loc=1
   if(is.na(end.loc)){
      truncated.chromosome.name <- as.character(sub("chr", "", chromosome))
      end.loc <- org.Hs.egCHRLENGTHS[[truncated.chromosome.name]]
      }

   if(!quiet){
      printf("requesting chipseq hits for%s in %s:%d-%d", tf, chromosome, start.loc, end.loc)
      }

   query <- sprintf("select * from chipseq where tf='%s' AND chrom='%s' AND start >= %d AND endpos <= %d",
                    tf, chromosome, start.loc, end.loc)
   tbl.tf <- dbGetQuery(db, query)

   dbDisconnect(db)

   if(!quiet){
      printf("chipseq database had %d hits in %d bases", nrow(tbl.tf), 1 + end.loc - start.loc)
      }

   tbl.tf.adjusted <- tbl.tf[, c("chrom", "peakStart", "peakEnd")]
   colnames(tbl.tf.adjusted) <- c("chrom", "start", "end")
   tbl.tf.adjusted$start <- tbl.tf[,"start"]
   tbl.tf.adjusted$end <- tbl.tf[,"endpos"]
   dim(tbl.tf.adjusted)


   tbl.motifs <- unique(requestMatchForRegions(fc, tbl.tf.adjusted, "hg38", fimoThreshold))
   if(!quiet){
      printf(" found %d fimo hits among %d ChIP-seq hits", nrow(tbl.motifs), nrow(tbl.tf.adjusted))
      }

   gr.chip <- GRanges(tbl.tf.adjusted)
   gr.motifs <- GRanges(tbl.motifs)
   tbl.overlaps <- as.data.frame(findOverlaps(gr.motifs, gr.chip, type="any"))
   colnames(tbl.overlaps) <- c("motif", "chip")
   chips.matched <- length(unique(tbl.overlaps$chip))
   chips.unmatched <- nrow(tbl.tf.adjusted) - chips.matched
   motifs.matched <- length(unique(tbl.overlaps$motif))
   fimo.pValScores <- -log10(tbl.motifs$pValue)
   dist <- fivenum(fimo.pValScores)

   if(!quiet){
      printf("found motif match for %d/%d ChIP hits (%5.2f%%)", chips.matched, nrow(tbl.tf.adjusted), 100 * chips.matched/nrow(tbl.tf.adjusted))
      printf("-log10(fimo pValue) range of %f %f %f %f %f",dist[1], dist[2], dist[3], dist[4], dist[5])
      }

  fimo.pValScores.withFailures <- c(fimo.pValScores, rep(-1, chips.unmatched))

  if(displayHistogram){
     hist(fimo.pValScores.withFailures,
     main=sprintf("JASPAR %s FIMO motif match to ChIP hits\n with %d failures (chr19)",tf, chips.unmatched),
     xlab="-log10(FIMO pValue)")
     } # if display

  invisible(tbl.motifs)

} # retreiveMatches
#----------------------------------------------------------------------------------------------------
test_CTCF_2kbRegionOnChr19 <- function()
{
  printf("--- test_CTCF_2kbRegionOnChr19")

  start <- 1590000
  end   <- 1592000
  tbl.motifs <- retrieveMatches("CTCF", "MA0139.1", "chr19", start.loc=start, end.loc=end,
                                displayHistogram=FALSE, fimoThreshold=1e-4, quiet=FALSE)
  checkEquals(dim(tbl.motifs), c(1,9))

} # test_CTCF_2kbRegionOnChr19
#----------------------------------------------------------------------------------------------------

