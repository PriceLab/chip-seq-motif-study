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
   bamFile <- "GSM749704_hg19_wgEncodeUwTfbsGm12878CtcfStdAlnRep1.bam"

   tbls <- identifyPeaks(bamFile, chrom="chr19", start=10000000, end=51190311)
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
#------------------------------------------------------------------------------------------------------------------------

