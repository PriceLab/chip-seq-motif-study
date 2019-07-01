library(FimoClient)
FIMO_HOST <- "khaleesi"
FIMO_PORT <- 60000

fc <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=FALSE)
start.loc <- 128486353
end.loc <- start.loc + 1000
tbl.regions <- data.frame(chrom="chr3", start=start.loc, end=end.loc, stringsAsFactors=FALSE)
tbl.matches <- requestMatchForRegions(fc, tbl.regions, "hg38", 1e-6)

head(tbl.matches)


