#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)

t=read.delim(paste0(prefix,".rates.tsv"),check.names=F)
t$Rate=factor(t$Rate,levels=as.character(t$Rate))
t=melt(t,id.vars="Rate")


for (rate in unique(t$Rate)) {
	png(paste0(prefix,".",rate,".rates.png"),width=length(unique(t$variable))*20+150)
	print(ggplot(t[t$Rate==rate,],aes(variable,value,color=variable))+
		geom_point()+
		theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
		coord_cartesian(ylim=c(0,max(t[t$Rate==rate,"value"])))+
		scale_color_discrete(guide=F)+
		xlab(NULL)+ylab("Rate"))
	dev.off()
}


t=read.delim(paste0(prefix,".ntrstat.tsv"),check.names=F,stringsAsFactors=F)
t$Condition=factor(t$Condition,levels=unique(t$Condition))

df=data.frame(
	Mode=rep(t$Mode,3),
	Type=rep(t$Type,3),
	Condition=rep(t$Condition,3),
	Parameter=rep(c("T->C","p_new","ntr"),each=nrow(t)),
	lower=c(t$`T->C`,t$p_new,t$ntr_lower),
	value=c(t$`T->C`,t$p_new,t$ntr),
	upper=c(t$`T->C`,t$p_new,t$ntr_upper)
)
df=df[!is.na(df$value),]

png(paste0(prefix,".ntr.png"),width=length(unique(t$Condition))*60+150)

ggplot(df,aes(Condition,value,color=Type,size=cut(upper-lower,breaks=c(0,0.01,0.05,0.1,Inf),include.lowest=TRUE)))+
	geom_point()+
	scale_size_manual("CI size",values=c(2,1,0.5,0.2,0.1))+
	facet_wrap(~Parameter,scales="free_y")+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
	
