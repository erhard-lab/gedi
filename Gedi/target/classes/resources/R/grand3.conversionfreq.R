#!/usr/bin/env Rscript

suppressPackageStartupMessages({
library(grandR)
library(ggplot2)
library(plyr)
})

	data=list(prefix=prefix)

	tab=read.delim(paste0(prefix,".conversion.freq.tsv.gz"))
	cond=unique(tab$Condition)

	
	subr=read.delim(paste0(prefix,".subread.tsv"))
	subr$Semantic=factor(as.character(subr$Semantic),levels=subr$Semantic)
	tab=merge(tab,subr,by="Subread")
	tab$Condition=factor(tab$Condition,levels=cond)
	ncond=nrow(tab)/nrow(unique(data.frame(tab$Category,tab$Genomic,tab$Read,tab$Semantic)))
	ncond=length(unique(tab$Condition))
	
	pdf(paste0(prefix,".conversion.freq.pdf"),width=if(ncond>=120) 10 else 5+ncond*1,height=7)
	for (cat in unique(tab$Category)) {
		print(PlotConversionFreq(data,category=cat)$plot)
	}
	g=dev.off()
	

	tab=read.delim(paste0(prefix,".conversion.knmatrix.tsv.gz"))
	cond=unique(tab$Condition)
	tab=merge(tab,subr,by="Subread")
	tab$Condition=factor(tab$Condition,levels=cond)
	ncond=length(unique(tab$Condition))
	tab$Covered=tab$n*tab$Count
	tab$Converted=tab$k*tab$Count
	tab=ddply(tab,.(Condition,Label,Semantic),function(s) data.frame(Covered=sum(s$Covered),Converted=sum(s$Converted)))
	
	pdf(paste0(prefix,".model.subreads.pdf"),width=3+ncond*0.5,height=7)
	for (lab in unique(tab$Label)) {
		t=tab[tab$Label==lab,]
		print(ggplot(t,aes(Condition,Covered,fill=Semantic))+
			cowplot::theme_cowplot()+
			geom_bar(stat="Identity")+
			scale_fill_brewer(NULL,palette="Dark2")+
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
		print(ggplot(t,aes(Condition,Covered,fill=Semantic))+
			cowplot::theme_cowplot()+
			geom_bar(stat="Identity",position='fill')+
			scale_fill_brewer(NULL,palette="Dark2")+
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
		print(ggplot(t,aes(Condition,Converted,fill=Semantic))+
			cowplot::theme_cowplot()+
			geom_bar(stat="Identity")+
			scale_fill_brewer(NULL,palette="Dark2")+
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
		print(ggplot(t,aes(Condition,Converted,fill=Semantic))+
			cowplot::theme_cowplot()+
			geom_bar(stat="Identity",position='fill')+
			scale_fill_brewer(NULL,palette="Dark2")+
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
	}
	g=dev.off()
	
