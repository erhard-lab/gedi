#!/usr/bin/env Rscript

suppressPackageStartupMessages({
	library(uniReg)
	library(ggplot2)
	library(reshape2)
	library(plyr)
	library(cowplot)
})

ks.test2=function(STATISTIC,n) exp(- 2 * n * STATISTIC^2)

lse=function(u,v) {
	m = max(u,v);
	log(exp(u-m)+exp(v-m))+m
}
decumsum=function(x) x-c(0,x[-length(x)])
rm.r0=function(m) m[rowSums(m)>0,]


file = paste0(prefix,".pep.fdrdata.tsv")
out= paste0(prefix,".pep.fdr.fit.tsv")
out2 = paste0(prefix,".pep.fdr.estimates.tsv")
out3 = paste0(prefix,".pep.fdr.test.tsv")


t=read.delim(file)

if (!("Genome" %in% names(t))) t$Genome="";


pdf(paste0(file,".pdf"))

re=data.frame(Length=c(),Annotation=c(),Genome=c(),ALC=c(),FDR=c(),PEP=c(),TargetDist=c())
re2=data.frame(Length=c(),Annotation=c(),Genome=c(),TrueTargets.lower=c(),TrueTargets=c(),TrueTargets.upper=c(),FalseTargets.lower=c(),FalseTargets=c(),FalseTargets.upper=c())


for (l in sort(unique(t$Length))) {
	
	tt=t[t$Length==l & t$Annotation=="Total",]
	
	if (sum(tt$Decoy)>=10) {
		cat(sprintf("Computing decoy model for l=%d...\n",l))
		fit=unireg(tt$ALC,tt$Decoy,sigmasq=1,nCores=nthreads)
		Fitted=pmax(0,fit$fitted.values)/sum(pmax(0,fit$fitted.values))
		CFitted=rev(cumsum(rev(pmax(0,fit$fitted.values))))/sum(pmax(0,fit$fitted.values))
	
		print(ggplot(tt,aes(ALC,rev(cumsum(rev(Decoy)))/sum(Decoy)))+geom_point()+ggtitle(sprintf("Decoy score distribution (l=%d)",l))+ylab("Cumulative freq.")+geom_line(aes(y=CFitted)))
	
		cat(sprintf("Computing target model for l=%d\n",l))
		tc=t[t$Length==l & t$Annotation=="Reference",]
		fit=unireg(tc$ALC,pmax(0,tc$Target-tc$Decoy),sigmasq=1,nCores=nthreads)
		pFitted=pmax(0,fit$fitted.values)/sum(pmax(0,fit$fitted.values))
		pCFitted=rev(cumsum(rev(pmax(0,fit$fitted.values))))/sum(pmax(0,fit$fitted.values))
		
		print(ggplot(tc,aes(ALC,rev(cumsum(rev(pmax(0,Target-Decoy))))/sum(pmax(0,Target-Decoy))))+geom_point()+ggtitle(sprintf("Target-Decoy score distribution (l=%d)",l))+ylab("Cumulative freq.")+geom_line(aes(y=pCFitted)))
	}
	
	tag = unique(t[,c("Annotation","Genome")])
	
	for (r in 1:dim(tag)[1]) {
		cat=tag[r,"Annotation"]
		genome=tag[r,"Genome"]
		tc=t[t$Length==l & t$Annotation==cat & t$Genome==genome,]
		
		if (sum(tc$Target+tc$Both)>0) {
			cat(sprintf("Computing mixture model for l=%d %s...\n",l,cat))
			
			loglik=function(x) sum(lse(log(x)+log(pFitted+1E-6),log1p(-x)+log(Fitted+1E-6))*(tc$Target+tc$Both))
			mod=optimize(loglik,interval=0:1,maximum=T)$maximum
			if (loglik(mod)-qchisq(0.95,1)/2<loglik(0)) lower=0 else lower=uniroot(function(z) loglik(mod)-loglik(z)-qchisq(.95,1)/2,c(0,mod))$root
			if (loglik(mod)-qchisq(0.95,1)/2<loglik(1)) upper=1 else upper=uniroot(function(z) loglik(mod)-loglik(z)-qchisq(.95,1)/2,c(mod,1))$root
			
			plot(0:100,tc$Target+tc$Both,pch=16,xlab="ALC",ylab="Targets",main=sprintf("Target fit (%s, %s, l=%d)",cat,genome,l))
			lines(0:100,pFitted*mod*sum(tc$Target+tc$Both)+Fitted*(1-mod)*sum(tc$Target+tc$Both),col='red',lwd=5)
			lines(0:100,pFitted*lower*sum(tc$Target+tc$Both)+Fitted*(1-lower)*sum(tc$Target+tc$Both),col='red',lwd=3,lty=2)
			lines(0:100,pFitted*upper*sum(tc$Target+tc$Both)+Fitted*(1-upper)*sum(tc$Target+tc$Both),col='red',lwd=3,lty=2)
			
			est.true=pFitted*mod*sum(tc$Target+tc$Both)
			est.false=Fitted*(1-mod)*sum(tc$Target+tc$Both)
			est.true.lower=pFitted*lower*sum(tc$Target+tc$Both)
			est.false.upper=Fitted*(1-lower)*sum(tc$Target+tc$Both)
			est.true.upper=pFitted*upper*sum(tc$Target+tc$Both)
			est.false.lower=Fitted*(1-upper)*sum(tc$Target+tc$Both)
			
			FDR=cummin(rev(cumsum(rev(est.false)))/rev(cumsum(rev(est.true+est.false))))
			PEP=ifelse(est.false+est.true==0,1,est.false/(est.true+est.false))
			
			targetdist=(tc$Target+tc$Both)*(1-PEP)
			targetdist=cumsum(targetdist)
		
			cumtargets=rev(cumsum(rev(tc$Target+tc$Both)))
			cumdecoys=rev(cumsum(rev(tc$Decoy)))
			targetcdf=pCFitted
			decoycdf=CFitted
			
			re=rbind(re,data.frame(Length=l,Annotation=cat,Genome=genome,ALC=0:100,FDR=FDR,PEP=PEP,TargetDist=targetdist,Cumulative.targets=cumtargets,Cumulative.decoys=cumdecoys,Target.cdf=targetcdf,Decoy.cdf=decoycdf))
			re2=rbind(re2,data.frame(Length=l,Annotation=cat,Genome=genome,TrueTargets.lower=sum(est.true.lower),TrueTargets=sum(est.true),TrueTargets.upper=sum(est.true.upper),FalseTargets.lower=sum(est.false.lower),FalseTargets=sum(est.false),FalseTargets.upper=sum(est.false.upper)))
		}	
	}
}


write.table(re,out,row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
write.table(re2,out2,row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)

getcutoffrow=function(tt,fdr=0.01) min(Inf,which(tt$FDR<fdr))

testtab=ddply(re,.(Annotation,Genome,Length),function(t) {
				cutoff=getcutoffrow(t); 
				if (is.finite(cutoff)) k=t$Cumulative.targets[cutoff] else k=0
				n=max(0,t$Cumulative.targets)
				if (is.finite(cutoff)) p=t$Decoy.cdf[cutoff] else p=0
				if (is.finite(cutoff)) lp=pbinom(k,size=n,prob=p,lower.tail=FALSE,log=TRUE) else lp=NA
		data.frame(All.Peptides=n,Decoy.p=p,Cutoff.Peptides=k,logpval=lp,minPeptides=qbinom(0.99,size=n,prob=p))
	})
write.table(testtab,out3,row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)



for (genome in setdiff(unique(re$Genome),'-')) {
	for (l in sort(unique(re$Length))) {
		df=re[re$Length==l & re$Genome==genome,]
		df=ddply(df,.(Annotation),function(s) cbind(s,data.frame(TargetDistNorm=s$TargetDist/max(s$TargetDist))))
		
		if (dim(df)[1]>0) {
			pvals=setNames(sapply(levels(df$Annotation), function(cat) {
								m=rm.r0(cbind(decumsum(df$TargetDist[df$Annotation=="Reference"]),decumsum(df$TargetDist[df$Annotation==cat])))
								if (!is.matrix(m) || any(dim(m)<2)) 1 else chisq.test(m)$p.value;
							}),levels(df$Annotation))
			levels(df$Annotation)=sprintf("%s (p=%.2g)",names(pvals),pvals)
			print(ggplot(df,aes(ALC,TargetDistNorm,color=Annotation))+geom_line()+ggtitle(sprintf("True targets score distribution (%s,l=%d)",genome,l))+theme(legend.position=c(0,1),legend.justification=c(0,1)))
		}
		
	}
	
#	for (l in sort(unique(re2$Length))) {
#		df=re2[re2$Length==l,]
#		print(ggplot(df,aes(TrueTargets,FalseTargets,col=Annotation))+geom_point(size=5)+xlab("True peptides")+ylab("False peptides")+ggtitle(sprintf("Estimated number of peptides (l=%d)",l)))
#	}


	percdf=data.frame(Annotation=c(),`CDS.percentile`=c(),cfreq=c())
	
	for (cat in setdiff(unique(re$Annotation),c("Total","Reference"))) {
	
		perc=matrix(ncol=2,nrow=0)
		
		for (l in sort(unique(re$Length))) {
			df=rbind(re[re$Length==l & re$Genome==genome,],re[re$Length==l & re$Annotation=="Reference",])
			if (sum(df$TargetDist>0) & sum(df$Annotation=="Reference")>0) {
				df=ddply(df,.(Annotation),function(s) cbind(s,data.frame(TargetDistNorm=s$TargetDist/max(s$TargetDist))))
				if (sum(df$Annotation==cat)>0) perc=rbind(perc,cbind(df$TargetDistNorm[ df$Annotation=='Reference' ],decumsum(df$TargetDist[ df$Annotation==cat ])))
			}
		}
		
		if (dim(perc)[1]>0) {
			perc=setNames(aggregate(perc[,2],list(perc[,1]),sum),c("CDS.percentile","cfreq"))
			perc[,2]=cumsum(perc[,2])/sum(perc[,2])
			
			n=sum(re$TargetDist[ re$Annotation==cat & re$ALC==100 ])
			stat=max(apply(perc,1,function(v) v[2]-v[1]))
			pval=ks.test2(stat,n)
			
			percdf=rbind(percdf,cbind(data.frame(Annotation=sprintf("%s p=%.2g",cat,pval)),perc))
		}
	}
	
	print(ggplot(percdf,aes(`CDS.percentile`,cfreq,color=Annotation))+geom_abline(linetype=2)+geom_line(size=1)+theme(legend.position=c(1,0),legend.justification=c(1,0))+ylab("Cumulative freq."))

}

invisible(dev.off())

