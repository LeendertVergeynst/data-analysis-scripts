require(xlsx)
data<-read.xlsx("PCA alkanes.xlsx",1)
X<-data[,-c(1:5)]
labs<-data$N

col<-factor(data$temp)
levels(col)<-c("black","blue","red","green","yellow")

pca1<-princomp(X,cor=F)
biplot(pca1)
pca1$loadings


par(mfrow=c(1,2))
plot(pca1$scores[,1],pca1$scores[,2],col="white")
text(pca1$scores[,1],pca1$scores[,2],labs,col=as.character(col))

plot(pca1$scores[,1],pca1$scores[,3],col="white")
text(pca1$scores[,1],pca1$scores[,3],labs,col=as.character(col))

summary(pca1)

