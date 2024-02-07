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


t=read.csv(gzfile(file),check.names=F)
t$Decoy=factor(revalue(t$Decoy,c(D="Decoy",T="Target",B="Ambiguous")),levels=c("Target","Decoy","Ambiguous"))
t$Annotation=revalue(t$Annotation,c(UTR5="5'-UTR",UTR3="3'-UTR"))
biggest = names(sort(table(t$Annotation),decreasing=TRUE))[1]

plot.genome=function(tt,geno) {

	df=ddply(tt,.(Annotation,Decoy),function(s) data.frame(qval=seq(0,0.1,by=0.001),n=dim(s)[1],cdf=ecdf(s$Q)(seq(0,0.1,by=0.001))))
	df=restrict.df(df)
	print(ggplot(df,aes(qval,cdf*n,col=Annotation,linetype=Decoy))+geom_vline(xintercept=0.01,linetype=2)+geom_line(size=1)+scale_color_brewer(NULL,palette="Set3")+scale_linetype_discrete(NULL)+ylab("Number of distinct peptides")+xlab("FDR")+ggtitle(geno))
	
	df=ddply(tt[tt$Annotation!=biggest,],.(Annotation,Decoy),function(s) data.frame(qval=seq(0,0.1,by=0.001),n=dim(s)[1],cdf=ecdf(s$Q)(seq(0,0.1,by=0.001))))
	df=restrict.df(df)
	if (dim(df)[1]>0) print(ggplot(df,aes(qval,cdf*n,col=Annotation,linetype=Decoy))+geom_vline(xintercept=0.01,linetype=2)+geom_line(size=1)+scale_color_brewer(NULL,palette="Set3")+scale_linetype_discrete(NULL)+ylab("Number of distinct peptides")+xlab("FDR")+ggtitle(geno))
	
	
	#df=ddply(tt,.(Annotation,Decoy),function(s) data.frame(qval=seq(0,0.1,by=0.001),n=dim(s)[1],cdf=ecdf(s$wrongQ)(seq(0,0.1,by=0.001))))
	#df=restrict.df(df)
	#print(ggplot(df,aes(qval,cdf*n,col=Annotation,linetype=Decoy))+geom_vline(xintercept=0.01,linetype=2)+geom_line(size=1)+scale_color_brewer(NULL,palette="Set3")+scale_linetype_discrete(NULL)+ylab("Number of distinct peptides")+xlab("wrong FDR")+ggtitle(geno))
	
	#df=ddply(tt[tt$Annotation!=biggest,],.(Annotation,Decoy),function(s) data.frame(qval=seq(0,0.1,by=0.001),n=dim(s)[1],cdf=ecdf(s$wrongQ)(seq(0,0.1,by=0.001))))
	#df=restrict.df(df)
	#if (dim(df)[1]>0) print(ggplot(df,aes(qval,cdf*n,col=Annotation,linetype=Decoy))+geom_vline(xintercept=0.01,linetype=2)+geom_line(size=1)+scale_color_brewer(NULL,palette="Set3")+scale_linetype_discrete(NULL)+ylab("Number of distinct peptides")+xlab("wrong FDR")+ggtitle(geno))
}


pdf(paste0(file,".pdf"))
t$Decoy[t$Decoy=='Ambiguous']="Target"

if ("Genome" %in% names(t)) {
	genomes = unique(t$Genome)
	genomes = genomes[!grepl(";",genomes)]
	
	for (geno in genomes) plot.genome(t[t$Genome==geno,],geno)
} else {
	plot.genome(t,"");
}



invisible(dev.off())
