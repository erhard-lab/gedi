#!/usr/bin/env Rscript

suppressPackageStartupMessages({
	library(grandR)
})

data=list(prefix=prefix)
file=paste0(prefix,".mismatch.raw.position.tsv.gz")

tab=read.delim(file,check.names=FALSE)

pdf(gsub(".tsv.gz",".bycondition.pdf",file),width=12,height=7)
for (cond in names(tab)[8:dim(tab)[2]]) {
	print(PlotMismatchPositionForSample(data,sample=cond)$plot)
}
g=dev.off()

if (dim(tab)[2]-7<=27) {
	pdf(gsub(".tsv.gz",".bytype.pdf",file),width=length(unique(tab$Sense))*5+2,height=length(unique(tab$Category))*2+0.2)
	print(PlotMismatchPositionForType(data,genomic="T",read="C")$plot)	
	print(PlotMismatchPositionForType(data,genomic="A",read="G")$plot)
	print(PlotMismatchPositionForType(data,genomic="G",read="A")$plot)	
	print(PlotMismatchPositionForType(data,genomic="C",read="T")$plot)	
	
	print(PlotMismatchPositionForType(data,genomic="T",read="A")$plot)	
	print(PlotMismatchPositionForType(data,genomic="T",read="G")$plot)	
	print(PlotMismatchPositionForType(data,genomic="A",read="C")$plot)	
	print(PlotMismatchPositionForType(data,genomic="A",read="T")$plot)	
	print(PlotMismatchPositionForType(data,genomic="G",read="C")$plot)	
	print(PlotMismatchPositionForType(data,genomic="G",read="T")$plot)	
	print(PlotMismatchPositionForType(data,genomic="C",read="A")$plot)	
	print(PlotMismatchPositionForType(data,genomic="C",read="G")$plot)	
	g=dev.off()
}
	
	

	
