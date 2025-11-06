library(plyr)
library(reshape2)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(ggseqlogo)

ggplot <- function(...) ggplot2::ggplot(...) + scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2")

plotTimecourse <- function(fileIn, fileOut) {
  df = read.delim(fileIn)

  df = ddply(df, c("DataSet", "ORF"), function(x) {
    dff = x
    dff$ReadCount = (dff$ReadCount/max(dff$ReadCount)) * 100
    return(dff)
  })

  df$Addition = factor(df$Addition, levels=c("Normal", "GCV"))
  df$ReadCountType = factor(df$ReadCountType, levels=c("Total RNA", "New RNA"))
  df = df[which(df$ORF != "-"),]
  pdf(fileOut)
  print(ggplot(df, aes(x=Timepoint,y=ReadCount,shape=Addition,color=DataSet,linetype=ReadCountType,group=paste(TiSS,DataSet,Addition,ReadCountType))) +
          geom_point() + geom_line() + facet_wrap(~ORF) + theme(legend.position="none"))
  dev.off()
}

plotTimecourse(fileIn, fileOut)