findBindingSites<- function(motifName, chrom, start, end)
{
  library(trena)
  library(MotifDb)
  
  pfms <- as.list(query (query(MotifDb, motifName), "jaspar2018"))
  #browser()
  motifMatcher <- MotifMatcher(genomeName="hg19", as.list(pfms), quiet=TRUE)
  
  tbl.regions <- data.frame(chrom, start, end, stringsAsFactors=FALSE)
  tbl.hits <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions, pwmMatchMinimumAsPercentage=85)
  }
#--------------------------------------------------------------------------------------
runTests <- function()
{
  test_findBindingSites()
  
  } # runTests
#------------------------------------------------------------------------------------------------------------------------
test_findBindingSites <- function()
{
  message(sprintf("--- test_findBindingSites"))
  tbl.hits<- findBindingSites(motifName="MA0139.1", chrom="chr19", start=53160025, end=56170311)
  checkTrue(is.list(tbl.hits))
  checkTrue(is.data.frame(tbl.hits))
  checkTrue(nrow(tbl.hits) > 0)
    
  } # test_findBindingSites
