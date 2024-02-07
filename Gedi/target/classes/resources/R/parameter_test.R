comp.p=function(x) sum(apply(x,2,function(c) sum(c*as.integer(rownames(x))))) / sum(apply(x,1,function(c) sum(c*as.integer(colnames(x)))))
#exp.x=function(x,p,oc) {for (n in 1:dim(x)[2]) x[1:2,n]=sum(x[3:(dim(x)[1]-1),n]/sum(dbinom(2:(dim(x)[1]-1),n,p)))*dbinom(0:1,n,p); x[is.na(x)]=0; x}
exp.x=function(x,p,oc) {
	if (!all(as.integer(rownames(x))==0:(dim(x)[1]-1))) stop("Incomplete!")
	for (c in 1:dim(x)[2]) {
		occ = oc[,c]
		k=as.integer(rownames(x))
		n=as.integer(colnames(x))[c]
		x[!occ,c]=sum(x[occ,c])/sum(dbinom(k[occ],n,p))*dbinom(k[!occ],n,p)
	}
	x[is.na(x)]=0
	x
}

TO=rbinom(9E6,50,0.3)
#c=rbinom(length(TO),TO,0.02)
c=rbinom(length(TO),TO,ifelse(runif(length(TO))<0.25,0.02,3E-4))
x=as.matrix(table(c,TO))


estimate=function(x, errp=5E-4) {
	onlyconv=apply(sapply(1:dim(x)[2],function(c) dbinom(as.integer(rownames(x)),as.integer(colnames(x)[c]),errp)*sum(x[,c]))<0.01*x,2,cummax)==1
	onlyconv=onlyconv&x>0
	onlyconv=onlyconv&matrix(rep(apply(onlyconv,2,sum)>1,each=dim(onlyconv)[1]),nrow=dim(onlyconv)[1])
	if (sum(onlyconv)==0) stop("Cannot estimate!")
	
	init.x=x
	init.x[matrix(rep(apply(onlyconv,2,sum)==0,each=dim(x)[1]),nrow=dim(x)[1])]=0

	
	inter=c(0,1)
	while(TRUE) {
		p=sum(inter)/2
		x=exp.x(x,p,onlyconv)
		np=comp.p(x)
		if (np<p) inter=c(inter[1],p) else inter=c(p,inter[2])
		print(c(inter,p,np))
		if (inter[2]-inter[1]<1E-12) break;
	}
	conv=sum(inter)/2
	x=init.x-x
	x[x<0]=0
	print(x)
	
	err1=comp.p(x)
	
	onlyconv[1,]=TRUE
	onlyconv[-1,]=FALSE
	inter=c(0,1)
	while(TRUE) {
		p=sum(inter)/2
		x=exp.x(x,p,onlyconv)
		np=comp.p(x)
		if (np<p) inter=c(inter[1],p) else inter=c(p,inter[2])
		print(c(inter,p,np))
		if (inter[2]-inter[1]<1E-12) break;
	}
	err=sum(inter)/2
	
	list(err=err,err1=err1,conv=conv)
}

