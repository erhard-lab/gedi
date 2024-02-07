#!/usr/bin/env Rscript

suppressPackageStartupMessages({
	library(ggplot2)
	library(cowplot)
	library(plyr)
})

restrict.df=function(df) {
	ddf=ddply(df,.(Annotation),function(s) max(s$cdf*s$n))
	use=ddf$Annotation[ddf$V1>max(ddf$V1)*0.01]
	df[df$Annotation %in% use,]
}


plot.genome=function(tt,geno) {
		print(ggplot(tt[tt$Annotation %in% c("CDS","5'-UTR","3'-UTR","ncRNA","OffFrame","Decoy"),],aes(`netMHC % rank`,color=Annotation))+stat_ecdf(size=1)+coord_cartesian(xlim=c(0,10))+scale_color_manual(NULL,values=cols)+ylab("Cumulative frequency")+ggtitle(geno))
		if (sum(tt$Annotation %in% c("Intergenic","Intronic"))>0) print(ggplot(tt[tt$Annotation %in% c("CDS","Intergenic","Intronic","Decoy"),],aes(`netMHC % rank`,color=Annotation))+stat_ecdf(size=1)+coord_cartesian(xlim=c(0,10))+scale_color_manual(NULL,values=cols)+ylab("Cumulative frequency")+ggtitle(geno))
		if (sum(tt$Annotation %in% c("PeptideSpliced","Substitution","Frameshift"))>0) print(ggplot(tt[tt$Annotation %in% c("CDS","PeptideSpliced","Substitution","Frameshift","Decoy"),],aes(`netMHC % rank`,color=Annotation))+stat_ecdf(size=1)+coord_cartesian(xlim=c(0,10))+scale_color_manual(NULL,values=cols)+ylab("Cumulative frequency")+ggtitle(geno))
		if (sum(tt$Annotation %in% c("RNASEQ","Extra"))>0) print(ggplot(tt[tt$Annotation %in% c("CDS","RNASEQ","Extra","Decoy"),],aes(`netMHC % rank`,color=Annotation))+stat_ecdf(size=1)+coord_cartesian(xlim=c(0,10))+scale_color_manual(NULL,values=cols)+ylab("Cumulative frequency")+ggtitle(geno))
		if (sum(tt$Annotation %in% c("Unknown"))>0) print(ggplot(tt[tt$Annotation %in% c("CDS","Unknown","Decoy"),],aes(`netMHC % rank`,color=Annotation))+stat_ecdf(size=1)+coord_cartesian(xlim=c(0,10))+scale_color_manual(NULL,values=cols)+ylab("Cumulative frequency")+ggtitle(geno))
		
		df=as.data.frame(table(Annotation=tt$Annotation,`HLA allele`=tt$`HLA allele`,Decoy=tt$Decoy))
		print(ggplot(df,aes(Decoy,Freq,fill=HLA.allele))+geom_bar(stat='identity',position='dodge')+facet_wrap(~Annotation,scales='free_y')+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle("1% FDR")+ggtitle(geno))
		print(ggplot(df,aes(Decoy,Freq,fill=HLA.allele))+geom_bar(stat='identity',position='fill')+facet_wrap(~Annotation)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle("1% FDR")+ggtitle(geno))
}

t=read.csv(file,check.names=F)
if ("netMHC % rank" %in% names(t) && !all(is.na(t$`netMHC % rank`))) {
	if ("Decoy netMHC % rank" %in% names(t)) {
		d=t[,!(names(t) %in% c("HLA allele","netMHC % rank"))]
		t=t[,!(names(t) %in% c("Decoy HLA allele","Decoy netMHC % rank"))]
		names(d)=names(t)
		
		t=ddply(t,.(Fraction,Scan),function(s) s[which(s$`netMHC % rank`==max(s$`netMHC % rank`))[1],])
		d=ddply(d,.(Fraction,Scan),function(s) s[which(s$`netMHC % rank`==max(s$`netMHC % rank`))[1],])
		d$Decoy="D"
		t$Decoy="T"
		t=rbind(t,d)
		t$UniqueAnnotation="Unidentified"
	}
	t=t[t$Q<0.01 | t$Decoy=='D',]
	t$Annotation=revalue(t$Annotation,c(UTR5="5'-UTR",UTR3="3'-UTR"))
	t$Decoy[t$Decoy=='B']='T'
	t$Annotation=factor(ifelse(t$Decoy=='D',"Decoy",as.character(t$Annotation)),levels=c(levels(t$Annotation),"Decoy"))
	cols=c(Decoy='red',setNames(sample(RColorBrewer::brewer.pal(11,"Spectral")),setdiff(levels(t$Annotation),"Decoy")))
	
	pdf(paste0(file,".netMHC.pdf"))
	
	if ("Genome" %in% names(t)) {
		for (geno in unique(t$Genome)) plot.genome(t[t$Genome==geno,],geno);
	} else {
		plot.genome(t,"");
	}
	
	invisible(dev.off())

}
