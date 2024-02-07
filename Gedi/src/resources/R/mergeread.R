
library(ggplot2)
library(reshape2)

t=do.call("rbind",lapply(files,function(f) cbind(Condition=gsub("^.*/(.*?).reads.tsv$","\\1",f),read.delim(f))))
t$Category=factor(t$Category,levels=unique(t$Category))
t=dcast(Condition~Category,data=t,value.var="Count")

write.table(t,out,row.names=F,col.names=T,sep="\t",quote=F)


for  (col in names(t)[-1]) {
	png(gsub(".tsv",sprintf(".%s.abs.png",col),out))
	print(ggplot(t,aes_string("Condition",paste0("`",col,"`")))+geom_bar(stat="identity")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
	dev.off()
}

t=cbind(Condition=t$Condition,t[,-(1:2)]/t[,2])
for  (col in names(t)[-1]) {
	png(gsub(".tsv",sprintf(".%s.rel.png",col),out))
	print(ggplot(t,aes_string("Condition",paste0("`",col,"`")))+geom_bar(stat="identity")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
	dev.off()
}

