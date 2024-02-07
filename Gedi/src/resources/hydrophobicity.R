#!/usr/bin/env Rscript

suppressPackageStartupMessages({
	library(ggplot2)
	library(cowplot)
	library(plyr)
})


plot.genome=function(tt,geno) {

	plot.ann=function(ttt,alpha,cols) {
		hqs=quantile(ttt$Hydrophobicity,c(0.01,0.99))
		print(ggplot(ttt,aes(RT,`Hydrophobicity`,color=Annotation,alpha=Annotation))+geom_point(size=1)+scale_alpha_manual(NULL,values=alpha)+coord_cartesian(ylim=hqs)+scale_color_manual(NULL,values=cols)+xlab("Retention time")+ylab("HI")+ggtitle(geno)+facet_wrap(~Fraction)+ guides(color = guide_legend(override.aes = list(alpha = 1))))
	}
	tt=tt[ order(tt$Annotation),]
	
	cols=c(Decoy='red',CDS="#3288bd",other='black')
	alpha=c(CDS=0.1,Decoy=0.1,other=1)
	plot.ann(tt[tt$Annotation %in% c("CDS","Decoy"),],alpha,cols)
	for (anno in setdiff(unique(tt$Annotation),c("CDS","Decoy"))) {
		names(alpha)[3]=anno
		names(cols)[3]=anno
		plot.ann(tt[tt$Annotation %in% c("CDS","Decoy",anno),],alpha,cols)
	}
}

t=read.csv(file,check.names=F)
if ("Hydrophobicity" %in% names(t) && !all(is.na(t$`Hydrophobicity`))) {
	t=t[t$Q<0.01 | t$Decoy=='D',]
	t$Annotation=revalue(t$Annotation,c(UTR5="5'-UTR",UTR3="3'-UTR"))
	t$Decoy[t$Decoy=='B']='T'
	t$Annotation=factor(ifelse(t$Decoy=='D',"Decoy",as.character(t$Annotation)),levels=c("Decoy",levels(t$Annotation)))
	
	pdf(paste0(file,".hydro.pdf"),height=4,width=4*length(unique(t$Fraction)))
	
	if ("Genome" %in% names(t)) {
		for (geno in unique(t$Genome)) plot.genome(t[t$Genome==geno,],geno);
	} else {
		plot.genome(t,"");
	}
	
	invisible(dev.off())

}
