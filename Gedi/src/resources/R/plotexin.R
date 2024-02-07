#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)

t=read.delim(paste0(prefix,".ext.tsv"))

t=melt(t[,grep("Gene|Readcount",names(t))],id.var="Gene")
levels(t$variable)=gsub(".Readcount","",levels(t$variable))
t=dcast(variable~Gene,data=t)

topic=setdiff(unique(gsub("Exonic|Intronic","",names(t)[-1])),"All")
   

pdf(paste(prefix,".exonintron.pdf",sep=''))

for (to in topic) {
if (all(paste0(c("Exonic","Intronic"),to) %in% names(t))) {
	tt=t[,c("variable",paste0(c("Exonic","Intronic"),to))]
	names(tt)=c("Condition","Exonic","Intronic")
	print(ggplot(tt,aes(Condition,Exonic/Intronic))+geom_bar(stat="identity")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(to))
}
}

dev.off()

