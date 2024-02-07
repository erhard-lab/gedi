#!/usr/bin/env Rscript

library(ggplot2)
library(plyr)
library(cowplot)

t=read.delim(file)
p=read.delim(param)

em = p$Value[p$Name=="Mean"]
esd = p$Value[p$Name=="Sd"]


pdf(paste(file,".pdf",sep=''))
ggplot(t,aes(Distance))+stat_ecdf()+coord_cartesian(xlim=c(0,600))+ylab("Cumulative frequency")+xlab("Distance annotated termination site to clone 5'")+stat_function(fun=pnorm,color='blue',size=2,geom="line",args=list(mean=em,sd=esd),n=600)+annotate(geom="text",x=100,y=0.1,label="Fit",size=12,color='blue')+ggtitle("All data")
t=t[t$Distance<450,]
ggplot(t,aes(Distance))+stat_ecdf()+coord_cartesian(xlim=c(0,600))+ylab("Cumulative frequency")+xlab("Distance annotated termination site to clone 5'")+stat_function(fun=pnorm,color='blue',size=2,geom="line",args=list(mean=em,sd=esd),n=600)+annotate(geom="text",x=100,y=0.1,label="Fit",size=12,color='blue')+ggtitle("Truncated data")
dev.off()

