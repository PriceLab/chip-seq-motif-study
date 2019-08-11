#Narrow peak Experimental data obtained from GEO (GSM3439936), G583 CTCF IP, replicate 1   ID: 303439936 (number 5 on the GEO list on Github)

unzippedTables <- lapply (c("GSM3439936_G583_ctcf_ip_rep1_WRT_G583_ctcf_input_rep1_peaks.narrowPeak","GSM1010734_hg19_wgEncodeHaibTfbsMcf7CtcfcV0422111PkRep1.broadPeak"), gzfile)
tables <- lapply(unzippedTables,read.table)
lapply(tables,head)
tbl.np <- tables[[1]]
tbl.bp <- tables[[2]]

colnames(tbl.np) <- c("chrom", "start", "end", "name", "score","width")
colnames(tbl.bp) <- c("chrom", "start", "end", "name", "score", "width")
tbl.np$width <- 1 + tbl.np$end - tbl.np$start
tbl.bp$width <- 1 + tbl.bp$end - tbl.bp$start

fivenum(tbl.np$width)
fivenum(tbl.bp$width)

boxplot(tbl.np$width, tbl.bp$width, names=c("Narrow Peaks","Broad Peaks"))
title(main="Comparison of DNA Region Width by Peak Type", xlab="Type of Peak", ylab="Width of DNA Region")

fivenum(tbl.np$score)
fivenum(tbl.bp$score)

tbl.bestNarrow <- subset(tbl.np, score == max(tbl.np$score))
tbl.bestBroad <- subset(tbl.bp, score == max(tbl.bp$score))
tbl.bestNarrow
tbl.bestBroad
