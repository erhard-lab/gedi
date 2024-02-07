xx=seq(0,25,by=.1)

x=rep(c(0,0.5,1,3,6,12,24),each=3)

plot.gene=function(n) {
	a=as.numeric(t[t$Gene==n,seq(104,144,by=2)])
	b=as.numeric(t[t$Gene==n,seq(105,145,by=2)])
	y=as.numeric(t[t$Gene==n,seq(36,96,by=3)])
	u=as.numeric(t[t$Gene==n,seq(37,97,by=3)])
	l=as.numeric(t[t$Gene==n,seq(35,95,by=3)])
	
	fun=function(p) -sum((a-1)*log(p[1]*exp(-p[2]*x))+(b-1)*log(1-p[1]*exp(-p[2]*x)))
	param=optim(c(0.7,5),fun)$par
	
	plot(x,y,col=rep(rainbow(3),8),ylim=c(0,1),main=n) 
	points(c(x,x),c(l,u),col=rep(rainbow(3),8),pch='-')
	lines(xx,param[1]*exp(-param[2]*xx))
	fun
}
fun=plot.gene(t[t$mESC.0h.chase.A.Readcount==152,1][1])
param=optim(c(0.9,4),fun)
plot(xx,exp(sapply(log(2)/xx,function(d) -fun(c(param$par[1],d)))),type='l')
