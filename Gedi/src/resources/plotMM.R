#!/usr/bin/env Rscript
# library(plyr)
# library(reshape2)
# library(ggplot2)

# cleanData <- function(df) {
#     tab = table(df[,3])
#     tab = tab[which(tab>=100)]
#     df = df[-(which(df[,3] %in% names(tab))),]
#     df = df[!grepl("K", df[,1]),]
#     df = df[!grepl("G", df[,1]),]
#     df = df[!grepl("M", df[,1]),]
#     return(df)
# }
#
# df = read.delim(dataFile)
# df <- cleanData(df)
# newfile = df[df[,3] > customThresh,]
# write.table(newfile, file=paste(prefix, "/manualSelection.tsv", sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
# df$x = 1:length(df[,1])
# df$y = sort(df[,3])
# pdf(pdfFile)
# print(ggplot(df, aes(x=x, y=y)) + geom_point() + geom_hline(aes(yintercept=threshMl, color="red")))
# dev.off()
