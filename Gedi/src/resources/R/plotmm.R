#!/usr/bin/env Rscript

library(ggplot2)
library(plyr)
library(data.table)

t=as.data.frame(fread(paste0(prefix,".mismatches.tsv")))
t$Mismatch=paste0(t$Genomic,"->",t$Read)
t$Mismatch=factor(t$Mismatch,levels=unique(t$Mismatch))
t$Rate=t$Mismatches/t$Coverage
t$se=sqrt(t$Rate*(1-t$Rate)/t$Coverage)
t$Condition=factor(t$Condition,as.character(unique(t$Condition)))
ncond=length(unique(t$Condition))

pdf(paste(prefix,".mismatches.pdf",sep=''),width=if(ncond>200) 7 else 4+length(unique(t$Condition)),height=7)

for (category in unique(t$Category)) {
	for (or in unique(t$Orientation)) {
		if (sum(t$Category==category & t$Orientation==or)>0) {
			if (ncond>200) {
				print(ggplot(t[t$Category==category & t$Orientation==or,],aes(Mismatch,Rate))+geom_boxplot()+theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste(category,or)))
			} else {
				print(ggplot(t[t$Category==category & t$Orientation==or,],aes(Mismatch,Rate,color=Condition))+geom_point(position=position_dodge(width=0.7))+geom_errorbar(aes(ymin=Rate-se,ymax=Rate+se),alpha=0.3,position=position_dodge(width=0.7),width=0)+theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste(category,or)))
			}
		}
	}
}
for (mm in unique(t$Mismatch)) {
	for (or in unique(t$Orientation)) {
		if (sum(t$Mismatch==mm & t$Orientation==or)>0) {
			if (ncond>200) {
				print(ggplot(t[t$Mismatch==mm & t$Orientation==or,],aes(Category,Rate))+geom_boxplot()+theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste(mm,or)))
			} else {
				print(ggplot(t[t$Mismatch==mm & t$Orientation==or,],aes(Category,Rate,color=Condition))+geom_point(position=position_dodge(width=0.7))+geom_errorbar(aes(ymin=Rate-se,ymax=Rate+se),alpha=0.3,position=position_dodge(width=0.7),width=0)+theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste(mm,or)))
			}
		}
	}
}	
dev.off()


t=as.data.frame(fread(paste0(prefix,".mismatchdetails.tsv")))
t$Mismatch=paste0(t$Genomic,"->",t$Read)
t$Mismatch=factor(t$Mismatch,levels=unique(t$Mismatch))
t$Rate=t$Mismatches/t$Coverage
#
#sp=ddply(t,.(Category,Condition,Position,Coverage,Genomic,Overlap),function(x) c(Coverage=sum(x$Coverage)))
## make s complete such that facet_wrap can be used
#
#pdf(paste(prefix,".readcoverage.pdf",sep=''),width=2+ceiling(sqrt(length(unique(t$Condition))))*5,height=2+10*ceiling(sqrt(length(unique(t$Condition)))))
#s=merge(sp[sp$Overlap==0,],data.frame(Category=rep(levels(sp$Category),length(levels(sp$Condition))),Condition=rep(levels(sp$Condition),each=length(levels(sp$Category)))),all=T)
#comp=is.na(s$Position)
#s[comp,"Position"]=0
#s[comp,"Genomic"]="A"
#s[comp,"Coverage"]=0
#s[comp,"Overlap"]=0
#print(ggplot(s,aes(Position,Coverage,color=Genomic))+geom_line()+facet_wrap(~paste(Condition,Category),ncol=length(unique(s$Category)),scales="free_y")+theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle("Non-overlap part"))
#s=merge(sp[sp$Overlap==1,],data.frame(Category=rep(levels(sp$Category),length(levels(sp$Condition))),Condition=rep(levels(sp$Condition),each=length(levels(sp$Category)))),all=T)
#comp=is.na(s$Position)
#s[comp,"Position"]=0
#s[comp,"Genomic"]="A"
#s[comp,"Coverage"]=0
#s[comp,"Overlap"]=1
#print(ggplot(s,aes(Position,Coverage,color=Genomic))+geom_line()+facet_wrap(~paste(Condition,Category),ncol=length(unique(s$Category)),scales="free_y")+theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle("Overlap part"))
#dev.off()
#
t$Overlap=factor(t$Overlap)
levels(t$Overlap)=paste0("Overlap=",levels(t$Overlap))
t$Opposite=factor(t$Opposite)
levels(t$Opposite)=paste0("Opposite=",levels(t$Opposite))

pdf(paste(prefix,".mismatchpos.pdf",sep=''),width=10,height=6)
for (category in unique(t$Category)) {
	print(ggplot(t[t$Category==category,],aes(Position,Rate,color=Read))+geom_line()+facet_grid(Overlap~Genomic+Opposite)+theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(category))
}
dev.off()


pdf(paste(prefix,".mismatchposzoomed.pdf",sep=''),width=10,height=6)
for (category in unique(t$Category)) {
	print(ggplot(t[t$Category==category,],aes(Position,Rate,color=Read))+geom_line()+facet_grid(Overlap~Genomic+Opposite)+coord_cartesian(ylim=c(0,0.04))+theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(category))
}
dev.off()

t=as.data.frame(fread(paste0(prefix,".doublehit.tsv")))
t$Mismatch=paste0(t$Genomic,"->",t$Read)
t$Mismatch=factor(t$Mismatch,levels=unique(t$Mismatch))
t$Rate=t$Hits/t$Coverage
t$se=sqrt(t$Rate*(1-t$Rate)/t$Coverage)
t$Condition=factor(t$Condition,as.character(unique(t$Condition)))

if (dim(t)[1]>0) {
	pdf(paste(prefix,".double.pdf",sep=''),width=4+length(unique(t$Condition)),height=7)
	
	for (category in unique(t$Category)) {
		if (ncond>200) {
			print(ggplot(t[t$Category==category,],aes(Mismatch,Rate))+geom_boxplot()+theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(category))
		} else {
			print(ggplot(t[t$Category==category,],aes(Mismatch,Rate,color=Condition))+geom_point(position=position_dodge(width=0.7))+geom_errorbar(aes(ymin=Rate-se,ymax=Rate+se),alpha=0.3,position=position_dodge(width=0.7),width=0)+theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(category))
		}
	}
	for (mm in unique(t$Mismatch)) {
		if (ncond>200) {
			print(ggplot(t[t$Mismatch==mm,],aes(Category,Rate))+geom_boxplot()+theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(mm))
		} else {
			print(ggplot(t[t$Mismatch==mm,],aes(Category,Rate,color=Condition))+geom_point(position=position_dodge(width=0.7))+geom_errorbar(aes(ymin=Rate-se,ymax=Rate+se),alpha=0.3,position=position_dodge(width=0.7),width=0)+theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(mm))
		}
	}	
	dev.off()
}
