#!/usr/bin/env Rscript

suppressPackageStartupMessages({
library(ggplot2)
library(cowplot)
library(grandR)
})
theme_set(theme_cowplot())

cumsumpercond=function(v,c,rel=FALSE) {
	re=rep(0,length(v))
	for (cc in unique(c)) re[c==cc]=if(rel) cumsum(v[c==cc])/sum(v[c==cc]) else cumsum(v[c==cc])
	re
}

	t=read.delim(paste0(prefix,".reads.lengths.tsv"),check.names=FALSE)
	s=read.delim(paste0(prefix,".reads.subreads.tsv"),check.names=FALSE)
	
	pdf(paste0(prefix,".targets.subreads.pdf"))
	ggplot(t,aes(Length,cumsumpercond(Frequency,Condition),col=Condition))+
		geom_step()+
		xlab("Overlapping read length")+ylab("Cumulative frequency")+
		theme(legend.position=c(1,0),legend.justification=c(1,0))
	
	ggplot(t,aes(Length,cumsumpercond(Frequency,Condition,rel=TRUE),col=Condition))+
		geom_step()+
		xlab("Overlapping read length")+ylab("Cumulative frequency")+
		theme(legend.position=c(1,0),legend.justification=c(1,0))
	
	ggplot(s,aes(Condition,Frequency,fill=Subread))+
		geom_bar(stat="identity",position="fill")+
		xlab("Positions in subread")+ylab("Frequency")
	
	ggplot(s,aes(Condition,Frequency,fill=Subread))+
		geom_bar(stat="identity")+
		xlab("Positions in subread")+ylab("Frequency")
	
	g=dev.off()
	
	
