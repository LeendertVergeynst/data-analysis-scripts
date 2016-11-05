test<-ptw.my(ref=t(dchroms2[,59]),samp=t(dchroms2[,46]),init.coef=c(0,1,0,0),optim.crit='WCC',mode='forward',verbose=T,trwdth=100,dt.max=edge(ref=ref,60,200))
plot((test$warp.fun-c(1:17950))[1,])

edge<-function(ref,dt.max,slope.length=100){
	e1<-seq(1,dt.max,length.out=slope.length)
	e2<-seq(dt.max,1,length.out=slope.length)
	c(e1,rep(dt.max,length(ref)-length(e1)-length(e2)),e2)
}

ref<-t(dchroms2[,59])
samp<-t(dchroms2[,46])
init.coef<-c(0,1,0,0)
optim.crit='WCC'
mode='forward'
verbose=T
trwdth=100


ref<-rfrnc
samp<-samp[i, , drop = FALSE]
optim.crit<-optim.crit
init.coef<-init.coef[i,]