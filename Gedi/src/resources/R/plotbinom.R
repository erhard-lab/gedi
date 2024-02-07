#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)

t=read.delim(paste0(prefix,".binomEstimated.tsv"))
t$Condition=factor(as.character(t$Condition),levels=as.character(t$Condition))

png(paste(prefix,".binomNewProportion.png",sep=''),width=800,height=800)
ggplot(t,aes(Condition,p_new*100,ymin=p_new_lower*100,ymax=p_new_upper*100,col=Status))+geom_point(size=1)+geom_errorbar(width=0.2)+ylim(0,100)+ylab("Percentage of new RNA")+theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

t=t[!grepl("no4sU",t$Condition),]
r=melt(t[,c(1:3,7)])

png(paste(prefix,".binomFreq.png",sep=''),width=800,height=800)
ggplot(r,aes(Condition,value,col=variable,shape=Status))+geom_point(size=2)+scale_color_brewer("",palette="Dark2")+theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("T->C Incorporation frequency\ninto new RNA")+xlab(NULL)
dev.off()


if (file.exists(paste0(prefix,".binomOverlapEstimated.tsv"))) {
	t=read.delim(paste0(prefix,".binomOverlapEstimated.tsv"))
	t$Condition=factor(as.character(t$Condition),levels=as.character(t$Condition))
	
	png(paste(prefix,".binomOverlapNewProportion.png",sep=''),width=800,height=800)
	print(ggplot(t,aes(Condition,p_new,ymin=p_new_lower,ymax=p_new_upper,col=Status))+geom_point(size=1)+geom_errorbar(width=0.2)+ylim(0,1)+ylab("Proportion of new RNA")+theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
	dev.off()

	t=t[!grepl("no4sU",t$Condition),]
	r=melt(t[,c(1:3,7)])
	png(paste(prefix,".binomOverlapFreq.png",sep=''),width=800,height=800)
	print(ggplot(r,aes(Condition,value,col=variable,shape=Status))+geom_point(size=2)+scale_color_brewer("",palette="Dark2")+theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("T->C Incorporation frequency\ninto new RNA")+xlab(NULL))
	dev.off()
	
}
