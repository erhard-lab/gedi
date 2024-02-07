#!/usr/bin/env Rscript

suppressPackageStartupMessages({
library(ggplot2)
library(patchwork)
})

tab=read.delim(paste0(prefix,".stats.tsv"),stringsAsFactors=TRUE)
tab$Mismatch=interaction(tab$Genomic,tab$Read,sep=">",drop=TRUE)


for (mm in levels(tab$Mismatch)) {
	g1=ggplot(tab[tab$Mismatch==mm,],aes(Position,Original))+
		geom_point(shape=1)+
		geom_point(shape=16,aes(y=Retained))+
		ylab("# mismatches")+
		xlab(NULL)+
		ggtitle(mm)+
		cowplot::theme_cowplot()
	g2=ggplot(tab[tab$Mismatch==mm,],aes(Position,Retained/Original))+
		geom_point()+
		scale_y_continuous("Retained",labels=scales::percent_format(1))+
		cowplot::theme_cowplot()
		
	png(paste0(prefix,".stats.plots/",prefix,".phred.",mm,".png"),width=max(tab$Position)*0.07+1,height=6,res=300,units="in")
	print(g1/g2)
	dev.off()
}

