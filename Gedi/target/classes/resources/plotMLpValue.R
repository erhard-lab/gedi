library(ggplot2)

multimapThresh <- 99999;

cleanData <- function(df) {
    tab = table(df[,3])
    tab = tab[which(tab>=multimapThresh)]
    todelete = which(df[,3] %in% names(tab))
    if (length(todelete) > 0) {
        df = df[-(todelete),]
    }
    return(df)
}

df = read.delim(dataFile)
df <- cleanData(df)
newfile = df[df[,3] < customThresh,]
write.table(newfile, file=paste(prefix, "/manualSelection.tsv", sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
df$x = 1:length(df[,1])
df$y = sort(df[,3])
pdf(pdfFile)
print(ggplot(df, aes(x=x, y=y)) + geom_point() + geom_hline(aes(yintercept=threshMl, color="red")))
dev.off()
