#!/usr/bin/env Rscript

library(ggplot2)
library(plyr)

t=read.delim(paste0(prefix,".errorEstimated.tsv"))

png(paste(prefix,".errorRates.png",sep=''),width=800,height=800)
ggplot(t,aes(Doublehit,Mismatch,color=Type))+geom_point()+geom_abline()+geom_abline(data=ddply(t,.(Type),function(s) data.frame(est=median(s$Mismatch-s$Doublehit))),aes(slope=1,intercept=est,color=Type),linetype=2)+xlab("Double-Hit rate")+ylab("Mismatch rate")+theme_bw()+theme(text=element_text(size=24))
dev.off()

