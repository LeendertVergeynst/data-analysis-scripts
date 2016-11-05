library(ptw)
library(xlsx)
root<-"~/Documents/15-16 ARC/Research/2015 diesel-sea water enrichment experiment/hydrocarbon analysis/GC-FID/"
source("functions.R")
source('my ptw.R')
sample.info<-read.xlsx(paste0(root,"sample info.xlsx"),1)
sample.info$path<-paste0(root,sample.info$folder,"/",sample.info$file)

# set subset1
chrom1<-read.table(sample.info$path[1],col.names=c('t','s'))
subset1<-seq(1,dim(chrom1)[1],by=4)
t1<-chrom1$t[subset1]

# load chroms1
chroms1<-sapply(sample.info$path,function(x) read.table2(x,subset1,2),simplify=T,USE.NAMES=F)
save(chroms1,file="chroms1.Rdata")


# remove not interesting part of chrom = subset 2
subset2		<-	c(3500:21000)
plot(chroms1[subset2,1],type='l',lwd=0.5)
for(i in 1:dim(chroms1)[2]) lines(chroms1[subset2,i],lwd=0.5)
lines(chroms1[subset2,11],lwd=0.5,col='red')
t2			<-	t1[subset2]
chroms2		<-	chroms1[subset2,]

# scale chroms

plot(chroms2[,8]/max(chroms2[,8]),type='l',lwd=0.5)
for(i in 1:dim(chroms1)[2]) lines(chroms2[,i]/max(chroms2[,i]),lwd=0.5)
lines(chroms2[,59]/max(chroms2[,59]),lwd=0.5,col='red')

maxchroms2	<-	apply(chroms2,2,function(x) max(x[1000:dim(chroms2)[1]]))
chroms3		<-	scale(chroms2,center=F,scale=maxchroms2)
startchroms3<-	chroms3[1,]
dchroms1		<-	apply(chroms3,2,diff)
add			<-	matrix(rep(0,dim(dchroms1)[2]*100),ncol=dim(dchroms1)[2])
dchroms2		<-	rbind(add,dchroms1,add)

# warping
plot(dchroms2[,1],type='l',lwd=0.5)
for(i in 1:dim(dchroms2)[2]) lines(dchroms2[,i],lwd=0.5)
lines(dchroms2[,59],lwd=0.5,col='red')

ptw1<-ptw.my(ref=t(dchroms2[,59]),samp=t(dchroms2[,]),init.coef=c(0,1,0,0,0,0,0),optim.crit='WCC',mode='forward',verbose=T,trwdth=400,dt.max=edge(ref=t(dchroms2[,59]),60,100))
plot(ptw1,what="function")
plot(ptw1,type="simultaneous",xlim=c(0000,5000))
plot(ptw1,type="simultaneous",xlim=c(10000,11000))
plot(ptw1,type="simultaneous",xlim=c(13000,15000))
plot(ptw1,type="simultaneous",xlim=c(15500,16500))
plot(ptw1,type="simultaneous",xlim=c(16500,17700))
lines(ptw1$reference[1,],col='blue')

ptw1$warped.sample[is.na(ptw1$warped.sample)]<-0

save(ptw1,file="ptw1.Rdata")

ptw2<-ptw.my(ref=apply(ptw1$warped.sample,2,median),samp=ptw1$warped.sample,init.coef=extend.coef(ptw1$warp.coef,0),optim.crit='WCC',mode='forward',verbose=T,trwdth=300,dt.max=edge(ref=apply(ptw1$warped.sample,2,median),40,100))
save(ptw2,file="ptw2.Rdata")

plot(ptw2,what="function")
plot(ptw2,type="simultaneous",xlim=c(000,2000))
plot(ptw2,type="simultaneous",xlim=c(2000,5000))
plot(ptw2,type="simultaneous",xlim=c(7000,11000))
plot(ptw2,type="simultaneous",xlim=c(10200,10700))
plot(ptw2,type="simultaneous",xlim=c(10000,11000))
plot(ptw2,type="simultaneous",xlim=c(15500,16500))
plot(ptw2,type="simultaneous",xlim=c(16500,17700))
lines(ptw2$reference[1,],col='blue')



ptw2.wcc<-wcc.cross(ptw2,400)
ptw1.wcc<-wcc.cross(ptw1,400)
boxplot(c(ptw1.wcc),c(ptw2.wcc),c(ptw2.wcc-ptw1.wcc))

boxplot(c(ptw2.wcc-ptw1.wcc))
median(c(ptw2.wcc-ptw1.wcc),na.rm=T)


temp<-apply(ptw1$warped.sample,1,cumsum)
offset<-matrix(rep(startchroms3,dim(temp)[1]),byrow=T,ncol=dim(temp)[2])
chroms4<-(temp+offset)*matrix(rep(maxchroms2,dim(temp)[1]),byrow=T,ncol=dim(temp)[2])
save(chroms4,file="chroms4.Rdata")

plot(chroms4[,1]/max(chroms4[,1]),type='l',lwd=0.5,ylim=c(0,1),xlim=c(15000,20000))
for(i in 1:dim(chroms4)[2]) lines(chroms4[,i]/max(chroms4[,i]),lwd=0.5)



### peak detection & integration
peaks<-lapply(1:118,function(x) pdi(chroms4[,x],base.int="slope"))

a<-pdi(chroms4[,98],plot=T,xlim=c(17000,17500),ylim=c(0,300000),base.int="slope")

alkanes<-peaks[[53]]
select<-order(alkanes$peak.area,decreasing=T)[1:23]
a<-pdi(chroms4[,53],plot=T,base.int="slope",xlim=c(10000,12000))
points(alkanes[select,"peakmax"],chroms4[alkanes[select,"peakmax"],53])

alk.time<-alkanes[sort(select),"peakmax"]
alk.name<-c("n-C8","n-C9","n-C10","n-C11","n-C12","n-C13","n-C14","n-C15","n-C16","n-C17","pritane","n-C18","phytane","n-C19","n-C20","n-C21","n-C22","n-C23","n-C24","n-C25","n-C26","n-C27","n-C28")
alk<-data.frame(alk.name,alk.time)
save(alk,file="alk.Rdata")
rt<-data.frame()
PA<-data.frame()

for(i in 1:118){
	for(j in 1:23){
		t<-alk[j,2]
		t1<-which(abs(peaks[[i]]$peakmax-t)==min(abs(peaks[[i]]$peakmax-t)))
		rt[j,i]<-peaks[[i]]$peakmax[t1]
		PA[j,i]<-peaks[[i]]$peak.area[t1]
	}
}

row.names(rt)<-alk[,1]
row.names(PA)<-alk[,1]

save(rt,file="rt.Rdata")
save(PA,file="PA.Rdata")


alk[,2]<-apply(rt[,c(5,53,114)],1,mean)
boxplot(t(rt-alk[,2]))
abline(h=c(-20,-15,-10,0,10,15,20),col='red')
PA2<-PA
PA2[abs(rt-alk[,2])>20]<-NA
write.xlsx(cbind(alk,PA2),'alkanes.xlsx',row.names=F)

PA3<-t(PA2)
PAc<-PA3
IS0<-

PA3[which(sample.info$concentration==0.1),"phytane"]
phytane.groupmeans<-aggregate(PA3[,"phytane"],by=list(sample.info$type),FUN=mean)

PA3<-cbind(PA3,IS=PA3[,"phytane"]/phytane.groupmeans[match(sample.info$type,phytane.groupmeans[,1]),2])
alk[,"IDL"]<-c(0.2,rep(0.1,20),0.2,0.2)

for(d in 1:3){
	for(i in 1:23){
		std<-sample.info$day==d & !is.na(sample.info$concentration)
		x<-sample.info$concentration[std]
		s<-which(x>=alk$IDL[i])
		x<-x[s]	
		y<-PA3[std,i][s]/PA3[std,"IS"][s]
		x<-x[!is.na(y)]
		y<-y[!is.na(y)]
		lm<-lm(y~x,weights=1/x)
		lm$coefficients[1]
		PAc[sample.info$day==d,i]<-(PA3[sample.info$day==d,i]-lm$coefficients[1])/lm$coefficients[2]
		png(paste0("figs/",colnames(PA3)[i],"_",d,".png"))
		par(mfrow=c(2,1))
		par(mar=c(4,4,1,1)+0.1)
		plot(x,y,pch=20,ylab="peak area")
		lines((x),predict(lm))
		plot(x,y/predict(lm)-1,pch=20,ylab="error")
		abline(h=0)
		text(0.75,0.05,paste0("RMSE = ",round(sqrt(sum((y/predict(lm)-1)^2)/(length(y)-2)),3)))
		dev.off()
	}
}

### inj peak detection and integration
plot(chroms1[c(3200:3500),8],type='l')
for(i in 1:dim(chroms1)[2]) lines(chroms1[c(3200:3500),i],lwd=0.5)

a<-pdi(chroms1[c(3200:3500),5],plot=T,span=0.05)
inj.peaks<-lapply(1:118,function(x) pdi(chroms1[c(3200:3500),x],span=0.05))
PA.inj.peaks<-unlist(lapply(inj.peaks, function(x) max(x$peak.area)))

PAc<-cbind(PAc,PA.inj.peak=PA.inj.peaks)

write.xlsx(PAc,"calibrated alkanes.xlsx")

