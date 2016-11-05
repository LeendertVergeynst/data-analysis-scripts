read.table2<-function(x,subset,col){
	t<-read.table(x,col.names=c('t','s'))[subset,col]
	return(t)
}

extend.coef<-function(coef,n){
	if(n==0) return(coef)
	cbind(coef,matrix(0,ncol=n,nrow=dim(coef)[1]))
}


wcc.cross<-function(x,trwdth=400){
	n.coef <- dim(x$warp.coef)[2]
	nt <- ncol(x$warped.sample)
	ns <- nrow(x$warped.sample)
	trima<-upper.tri(matrix(1,nrow=ns,ncol=ns), diag=F)*1
	trima[trima==0]<-NA
	B <- matrix(1:nt, nrow = nt, ncol = n.coef)
	B <- t(apply(B, 1, cumprod))/B
	x$warped.sample[is.na(x$warped.sample)]<-0
	for(i in 1:ns){
		for(j in 1:ns){
			if(is.na(trima[i,j])) next
			trima[i,j]<-1-wcc(x$warped.sample[i,], x$warped.sample[j,], trwdth, wghts = 1 - (0:trwdth)/trwdth)
		}
	}
	return(trima)
}

argmax <- function(y,w=20,...) {
  require(zoo)
  x<-1:length(y)
  n<-length(y)
  y.smooth <- loess(y ~ x,...)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  return(list(i=i.max,y.smooth=y.smooth))
}

pdi<-function(y,plot=F,span=0.001,base.int="flat",...){ #peak detection and integration
	y2<-argmax(y,w=30,span=span)
	peakmax<-y2$i
	peakmin<-argmax(-y,w=2,span=span)$i
	peaks<-data.frame(peakmax=peakmax)
	all.false<-function(x){
		if(sum(x*1)==0) return(NA)
		x
	}
	peaks[,"start"]<-sapply(peakmax,function(x) max(all.false(peakmin[peakmin<c(x)])))
	peaks[,"stop"]<-sapply(peakmax,function(x) min(all.false(peakmin[peakmin>c(x)])))
	peaks<-na.omit(peaks)
	rownames(peaks)<-NULL
	if(plot) plot(y,type='l',...)
	if(plot) lines(peaks$peakmax,y[peaks$peakmax],type='h',col='red',lty=2,lwd=0.5)
	if(plot) lines(y2$y.smooth,type='l',col='green',lwd=0.5)
	for(i in 1:dim(peaks)[1]){
		x0<-peaks[i,c("start")]
		x1<-peaks[i,c("stop")]
		base<-min(y[c(x0,x1)])
		peaks[i,"base"]<-base
		profile<-y[x0:x1]
		#peaks[i,'peakmax']<-which(profile==max(profile))+x0-1
		if(base.int=="flat") peaks[i,"peak.area"]<-sum(profile)-length(profile)*base
		if(base.int=="slope") peaks[i,"peak.area"]<-sum(profile)-length(profile)*sum(y[c(x0,x1)])/2
		if(plot & base.int=="flat") polygon(c(x0,x0:x1,x1),c(base,profile,base), col=rgb(1,0,0,0.3),border=NA)
		if(plot & base.int=="slope") polygon(c(x0,x0:x1,x1),c(y[x0],profile,y[x1]), col=rgb(1,0,0,0.3),border=NA)
	}
	return(peaks)
}