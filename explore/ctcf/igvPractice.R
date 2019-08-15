library(igvR) #you are opening a package
igv <- igvR() #connects to browser
Sys.sleep(2)

setGenome(igv, "hg19")
showGenomicRegion(igv, "CTCF")
Sys.sleep(5) # give igv.js time to load

tbl.narrowp <- data.frame(chrom=rep("chr14", 3), #best scoring chIP-seq region
                          start=106873431,
                          end=106873939,
                          name=c("a", "b", "c"),
                          score=runif(3),
                          strand=rep("*", 3),
                          stringsAsFactors=FALSE)
tbl.broadp <- data.frame(chrom=rep("chr11", 3), #best scoring chIP-seq region
                         start= 189622,
                         end=190360,
                         name=c("a", "b", "c"),
                         score=runif(3),
                         strand=rep("*", 3),
                         stringsAsFactors=FALSE)
NarrowpTrack <- DataFrameAnnotationTrack("track.2", tbl.narrowp, color="red", displayMode="SQUISHED")
BroadpTrack <- DataFrameAnnotationTrack("track.3", tbl.broadp, color="green", displayMode="SQUISHED")
displayTrack(igv, NarrowpTrack)
displayTrack(igv, BroadpTrack)
