library(igvR)

igv <- igvR()
setGenome(igv, "hg19")
showGenomicRegion(igv, "PDCD1")
chromLoc <- getGenomicRegion(igv)
shoulder <- 10000
bigStart <- chromLoc$start - shoulder
bigEnd   <- chromLoc$end + shoulder

biggerLoc <- with(chromLoc, sprintf("%s:%d-%d", chrom, start-shoulder, end+shoulder))
showGenomicRegion(igv, biggerLoc)
tbl.bed <- read.table(file="ENCFF001VMZ-small.bed", sep="\t", as.is=TRUE, nrow=-1, header=TRUE, stringsAsFactors=TRUE)
colnames(tbl.bed) <- c("chrom", "start", "end", "strand", "score", "x1", "x2", "x3", "x4", "x5")
dim(tbl.bed)
head(tbl.bed)

tbl.hits <- subset(tbl.bed, chrom=="chr2" & start >= bigStart & end <= bigEnd)
track <- DataFrameQuantitativeTrack("CTCF", tbl.hits[, c(1:3, 5)], autoscale=TRUE, color="red")
displayTrack(igv, track)


library(trena)
pfms <- query(MotifDb, c("jaspar2018", "CTCF", "MA0139"))
mm <- MotifMatcher("hg19", as.list(pfms), quiet=TRUE)
tbl.match <- findMatchesByChromosomalRegion(mm, tbl.hits, pwmMatchMinimumAsPercentage=75)

track <- DataFrameQuantitativeTrack("MA0139", tbl.match[, c("chrom", "motifStart", "motifEnd", "motifRelativeScore")],
                                    color="blue", autoscale=TRUE)
displayTrack(igv, track)

tbl.all <- data.frame(chrom="chr2", start=min(tbl.hits$start), end=max(tbl.hits$end), stringsAsFactors=FALSE)

tbl.match.all <- findMatchesByChromosomalRegion(mm, tbl.all, pwmMatchMinimumAsPercentage=75)

track <- DataFrameQuantitativeTrack("MA0139.all", tbl.match.all[, c("chrom", "motifStart", "motifEnd", "motifRelativeScore")],
                                    color="magenta", autoscale=TRUE)
displayTrack(igv, track)
