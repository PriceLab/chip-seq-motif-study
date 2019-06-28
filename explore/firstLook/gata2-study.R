library(RPostgreSQL)
db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="hg38", host="khaleesi")
dbListTables(db) # [1] "gtf" "motifsgenes" "chipseq"

tbl.gata1 <- dbGetQuery(db, "select * from chipseq where tf='GATA1'")
dim(tbl.gata1)

tbl.gata1 <- dbGetQuery(db, "select * from chipseq where tf='GATA1'")
# the gata2 gene, thought to be repressed by gata1
query <- "select * from chipseq where tf='GATA1' AND chrom='chr3' AND start >= 128476675 AND endpos <= 128495937"
tbl.gata1.near.gata2 <- dbGetQuery(db, query)
dim(tbl.gata1.near.gata2)

#--------------------------------------------------------------------------------
# now we want to display these ChIP-seq hits in the igv genome browser
#--------------------------------------------------------------------------------
library(igvR)
igv <- igvR()
setGenome(igv, "hg38")
showGenomicRegion(igv, "GATA2")

tbl.track <- tbl.gata1.near.gata2[, c("chrom", "start", "endpos", "name")]
track <- DataFrameAnnotationTrack("GATA1.tfbs", tbl.track, color="maroon")
displayTrack(igv, track)

