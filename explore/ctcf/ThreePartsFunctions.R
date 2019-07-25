#fimoTest
test_FIMO <- function() 
  printf("--- test_FIMO")

  tf <- "CTCF"
  ID <- "MA0139.1"
  chrom <-"chr19"
  
  start <- 1590000
  end <- 1592000
  
  fimoThreshold<- 1e-4
  
  tbl.motifs <- retrieveMatches(tf, ID, chrom, start.loc=start, end.loc=end,
                              displayHistogram=FALSE, fimoThreshold, quiet=FALSE)
checkEquals(dim(tbl.motifs), c(1,9))