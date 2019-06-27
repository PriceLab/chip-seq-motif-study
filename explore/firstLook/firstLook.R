library(RPostgreSQL)
db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="chipatlas", host="khaleesi")
dbListTables(db) # [1] "experiments" "peaks"
dbGetQuery(db, "select * from experiments limit 4")
dbGetQuery(db, "select * from peaks limit 4")

tf <- "PARK7"
tf <- "ZNF263"
tf <- "LXRalpha"
tf <- "MEF2C"
tbl <- dbGetQuery(db, sprintf("select * from experiments where antigen = '%s'", tf))
dim(tbl)

id <- unique(tbl$id)
tbl.hits <- dbGetQuery(db, sprintf("select * from peaks where experiment='%s'", id))
dim(tbl.hits)


db2 <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="hg38", host="khaleesi")
dbGetQuery(db2, "select count(*) from chipseq")  # 80,143,591

dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="chipatlas", host="khaleesi")
dbListTables(db2) # "gtf"         "motifsgenes" "chipseq"
dbGetQuery(db2, "select * from chipseq limit 5")

tf <- "KLF1"
tf <- "PARK7"
tf <- "ZNF263"
tf <- "LXRalpha"
tf <- "MEF2C"
tbl <- dbGetQuery(db2, sprintf("select * from chipseq where tf='%s'", tf))
dim(tbl)

library(trena)
pfms <- query(MotifDb, c("jaspar2018", "MEF2C"))
mm <- MotifMatcher("hg38", as.list(pfms), quiet=TRUE)

tbl.fixed <- tbl[, c("chrom", "peakStart", "peakEnd")]
colnames(tbl.fixed) <- c("chrom", "start", "end")
tbl.fixed$start <- tbl.fixed$start - 10
tbl.fixed$end   <- tbl.fixed$end + 10
tbl.match <- findMatchesByChromosomalRegion(mm, head(tbl.fixed, n=100), pwmMatchMinimumAsPercentage=70)
dim(tbl.match)

tbl.fixedBigger <- tbl[, c("chrom", "start", "endpos")]
colnames(tbl.fixedBigger) <- c("chrom", "start", "end")
dim(tbl.fixedBigger)
tbl.match <- findMatchesByChromosomalRegion(mm, head(tbl.fixedBigger, n=100), pwmMatchMinimumAsPercentage=75)
dim(tbl.match)



# library(FimoClient)

