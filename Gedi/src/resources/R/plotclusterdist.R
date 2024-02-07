#!/usr/bin/env Rscript


suppressPackageStartupMessages({
	library(ggplot2)
	library(reshape2)
	library(plyr)
	library(cowplot)
})
theme_set(theme_cowplot())

t=read.delim(paste0(prefix,".3p.tsv"))
fitn=setNames(read.delim(paste0(prefix,".3p.parameter"))$Value,c("mean","sd"))


png(paste(prefix,".3p.distr.png",sep=''))
ggplot(t,aes(Distance,cumsum(Frequency)/sum(Frequency)))+geom_line()+stat_function(fun=pgamma,args=list(shape=fitn[1]^2/fitn[2]^2,scale=fitn[2]^2/fitn[1]),color='red')+stat_function(fun=pnorm,args=as.list(fitn),color='blue')+ylab("Cumulative frequence")
dev.off()
