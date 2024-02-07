#!/usr/bin/env Rscript

library(ggplot2)
library(cowplot)

rpk=t[,grepl("Readcount",names(t))]/(t$Length/1000)
names(rpk)=gsub("Readcount","RPK",names(rpk))
t=cbind(t,rpk)

scale=colSums(rpk)/1E6

gene.plot=function(t,name) {
	g=t[t$Gene==name,]
	rpk=as.numeric(g[grepl("RPK",names(t))])
	tpm=rpk/scale
	a=as.numeric(g[grepl("alpha",names(t))])
	b=as.numeric(g[grepl("beta",names(t))])
	n=gsub(".alpha","",names(g)[grepl("alpha",names(t))])
	n=factor(n,levels=n)
	
	df=data.frame(a=a,b=b,n=n,mean=a/(a+b),lower=qbeta(0.05,a,b),upper=qbeta(0.95,a,b),tpm=tpm)
	df2=data.frame(n=n,RNA=rep(c("new","old"),each=length(tpm)),
					tpm=c(tpm*a/(a+b),tpm*b/(a+b)),
					lower=c(qbeta(0.05,a,b)*tpm,qbeta(0.05,b,a)*tpm), 
					upper=c(qbeta(0.95,a,b)*tpm,qbeta(0.95,b,a)*tpm))
	
	p1=ggplot(df,aes(n,mean,ymin=lower,ymax=upper))+geom_errorbar(width=0.2)+geom_point()+theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("Proportion of new RNA")+xlab(NULL)
	#p2=ggplot(df,aes(n,tpm))+geom_bar(stat="identity")+theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("TPM")+xlab(NULL)
	#p3=ggplot(df,aes(n,tpm*mean,ymin=tpm*lower,ymax=tpm*upper))+geom_bar(stat="identity")+geom_errorbar(width=0.2)+theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("NTPM")+xlab(NULL)
	#p4=ggplot(df,aes(n,tpm*(1-mean),ymin=tpm*(1-upper),ymax=tpm*(1-lower)))+geom_bar(stat="identity")+geom_errorbar(width=0.2)+theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("OTPM")+xlab(NULL)
	#p2=ggplot(df,aes(n,tpm))+geom_bar(stat="identity",color='black',fill='white')+theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("TPM")+xlab(NULL)+geom_bar(aes(y=mean*tpm),stat="identity",color='black')+geom_errorbar(width=0.2,aes(ymin=tpm*lower,ymax=tpm*upper))
	
	p2=ggplot(df2,aes(n,tpm,ymin=lower,ymax=upper,color=RNA))+geom_errorbar(width=0.2,alpha=0.5)+geom_point(size=2)+theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position=c(0.99,0.01),legend.justification=c(1,0))+ylab("TPM")+xlab(NULL)+scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x)))+scale_color_grey(NULL,start=0,end=0.7)+annotation_logticks(sides = "l")
	
	
	g=plot_grid(p1,p2,ncol=2,labels=c("A","B"),label_size=28)
	title <- ggdraw() + draw_label(name, size=28)
	g=plot_grid(title, g, ncol=1, rel_heights=c(0.1, 1))
	
	#ggsave("name.png",height=19,width=13,plot=g)
	list(plot=g,data=df)
}

pdf("ercc.pdf",width=14)
for (i in grep("ERCC",t$Gene)) {
	print(gene.plot(t,t$Gene[i])$plot)
}
dev.off()
	

pdf("mcmv.pdf",width=14)
for (i in which(!grepl("ENSMUSG",t$Gene) & !grepl("ERCC",t$Gene))) {
	print(gene.plot(t,t$Gene[i])$plot)
}
dev.off()

	

	