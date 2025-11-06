library(ggplot2)

multimapThresh <- 99999;

# this will clean the data of values that occure more than multimapThresh times.
# on default settings none will be romved
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
df = cleanData(df)
df$x = 1:length(df[,1])
df$up = sort(df[,3])
df$down = sort(df[,4])

newfile = df[(df[,3] > customThreshUp & df[,4] > customThreshDown),]
write.table(newfile, file=paste(manualSelectionFile, sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
if (length(df[,1]) > 10000) {
    png(paste(manualSelectionPdf, ".png", sep=""), width=1024, height=768)
} else {
    pdf(manualSelectionPdf)
}
print(ggplot(df, aes(x=x, y=up)) + geom_point() + geom_hline(aes(yintercept=customThreshUp, color="red")))
print(ggplot(df, aes(x=x, y=down)) + geom_point() + geom_hline(aes(yintercept=customThreshDown, color="red")))
dev.off()
if (length(df[,1]) > 10000) {
    png(paste(pdfFile, ".png", sep=""), width=1024, height=768)
} else {
    pdf(pdfFile)
}
print(ggplot(df, aes(x=x, y=up)) + geom_point() + geom_hline(aes(yintercept=upThreshAutoParam, color="red")))
print(ggplot(df, aes(x=x, y=down)) + geom_point() + geom_hline(aes(yintercept=downThreshAutoParam, color="red")))
dev.off()
