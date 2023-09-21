## PCA and Fst by chomromsome 

library(data.table)
library(ggplot2)

a1f<-list.files(pattern="ad1_filt")
a2f<-a1f
a2f<-gsub("ad1","ad2",a2f)
N<-length(a1f)
ids<-read.table("IDs.txt",header=FALSE)
temp<-gsub("ad1_filt_lycSpecPool_chrom","",a1f)
chrom<-gsub(".txt","",temp)

## Fst
P<-vector("list",24)
n<-vector("list",24) 
H<-vector("list",24)
for(i in 1:N){
	a1<-as.matrix(fread(a1f[i],header=F))
	a2<-as.matrix(fread(a2f[i],header=F))
	n[[i]]<-a1+a2
	P[[i]]<-a2/(a1+a2) ## non-ref
	p<-a2/(a1+a2) ## non-ref
	p[is.na(p)]<-0.001
	H[[i]]<-2*p*(1-p)
}
	
save(list=ls(),file="fst.rdat")

## compute all of the pairwise Fst by pair and chromosome (not for each SNP)
Npop<-27
Nx<-(Npop*(Npop-1))/2
fstGw<-matrix(NA,nrow=Nx,ncol=24)
fst90<-matrix(NA,nrow=Nx,ncol=24)
x<-1
for(i in 1:(Npop-1)){
	for(j in (i+1):Npop){
		for(k in 1:24){
			keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
			pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
			Ht<-2*pbar*(1-pbar)
			Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
			
			fstGw[x,k]<-mean(Ht-Hs)/mean(Ht)
			fst<-(Ht-Hs)/Ht
			fst90[x,k]<-quantile(fst,.9,na.rm=TRUE)
		}
		x<-x+1
	}
}

head(fst)

###############################################################
mord<-c(1,12,18:24,2:11,13:17)
fstGw.new<-fstGw[,mord]
fstauto<-data.frame(matrix(nrow=351,ncol=2))

x=1

for(i in 1:(Npop-1)){
  for(j in (i+1):Npop){
    
    fstauto[x,1]<-paste(ids[i,1]," x ",ids[j,1],sep="")
    fstauto[x,2]<-mean(fstGw.new[x,-24])
    
    x=x+1
    
  }
}

head(fstauto)

hist(fstauto$X2)

p<-seq(0,1,0.1)
q<-seq(0,1,0.1)
b<-lm(p~q)

plot(fstauto$X2,fstGw.new[,24])+
  abline(b)

fstGw.new2<-data.frame(fstauto,X=fstGw.new[,24])
colnames(fstGw.new2)<-c("pair","auto","Z")
head(fstGw.new2)

ggplot(data=fstGw.new2,aes(x=auto,y=Z))+
  xlab("Fst on autochromosome")+ylab("Fst on Z chromosome")+
  geom_point()+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)+
  geom_abline(slope = coef(b)[["q"]], 
              intercept = coef(b)[["(Intercept)"]])+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=24))


##### relationship between migration and the degree of Fst difference between X and auto ######
######### calculating migration Nm ########
hist(fstGw.new2$Z-fstGw.new2$auto)
fstGw.new2$Nm<-((1/fstGw.new2$auto)-1)/4
ggplot(data=fstGw.new2,aes(x=Nm,y=(Z-auto)/auto))+
  xlab("Degree of Migration: Nm")+ylab("Fst (Z-auto)")+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=24))

cor.test(fstGw.new2$Nm,(fstGw.new2$Z-fstGw.new2$auto)/fstGw.new2$auto)

ggplot(data=fstGw.new2,aes(x=auto,y=(Z-auto)))+
  xlab("Fst auto")+ylab("Fst (Z-auto)")+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=24))

cor.test(fstGw.new2$auto,(fstGw.new2$Z-fstGw.new2$auto)/fstGw.new2$auto)
### relationship between chromosome length and  Fst ############

chrome.length<-c()

for(i in 1:24){
  chrome.length[i]<-length(n[[i]][,1])
}
chrome.length<-chrome.length[mord]

pdf("Chromelength.cor.pdf",width=16,height=16)
par(mfrow=c(4,3))
par(mar=c(4.5,5.5,2.5,1.5))

for(i in 1: 12){

plot(chrome.length[-24],fstGw.new[i,-24],ylab="Fst",xlab="chrome length",cex.lab=2)
text(200000,quantile(fstGw.new[i,-24],0.8),ids[i,1],cex=2)

}

dev.off()
mord
##############################################################################
fstC<-c()
pairsd<-c()
pairsv<-c()
x<-1
for(i in 1:(Npop-1)){
  for(j in (i+1):Npop){
    fstC<-c()
    for(k in c(1:16,18:24)){
      Nloci<-length(P[[k]][,1])
      keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
      pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
      Ht<-2*pbar*(1-pbar)
      Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
      fstC<-c(fstC,(Ht-Hs)/Ht)
    }
    pairsd[x]<-sd(fstC,na.rm=TRUE)
    pairsv[x]<-pairsd[x]/sqrt(mean(fstC,na.rm=TRUE))
    x=x+1
  }
}

plot(fstGw.new2$auto,pairsv)
plot(fstGw.new2$auto,pairsv,pch=16,
     xlab="Fst across autochromosome", ylab="Standard error of Fst across autochromosome",
     cex.lab=1.5)

#############################################################################


fstC<-c()
pairsd<-c()
for(i in 1:(Npop-1)){
  for(j in (i+1):Npop){
    fstC<-c()
    for(k in 1:24){
  Nloci<-length(P[[k]][,1])
  keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
  pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
  Ht<-2*pbar*(1-pbar)
  Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
  fstC<-c(fstC,(Ht-Hs)/Ht)
    }
    pairsd<-sd(fstC)
  }
}


fst90<-matrix(NA,nrow=Nx,ncol=24)
x<-1
for(i in 1:(Npop-1)){
  for(j in (i+1):Npop){
    for(k in 1:24){
      keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
      pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
      Ht<-2*pbar*(1-pbar)
      Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
      
      fstGw[x,k]<-mean(Ht-Hs)/mean(Ht)
      fst<-(Ht-Hs)/Ht
      fst90[x,k]<-quantile(fst,.9,na.rm=TRUE)
    }
    x<-x+1
  }
}

####################################################################
## plot of mean and 90 quantile per chromsome for each pair
mord<-c(1,12,18:24,2:11,13:17)
pdf("FstLycSpc.pdf",width=8,height=10)
par(mfrow=c(5,3))
par(mar=c(4,5,2.5,.5))
x<-1
for(i in 1:(Npop-1)){
	for(j in (i+1):Npop){
		plot(fstGw[x,mord],pch=19,ylim=c(0,1),xlab="Chromosome",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=.9)
		segments(x0=1:24,y0=fstGw[x,mord],x1=1:24,y1=fst90[x,mord])
		title(main=paste(ids[i,1]," x ",ids[j,1],sep=""),cex.main=1.2)
		mn<-round(mean(fstGw[x,mord]),2)
		text(5,.9,mn)
		x<-x+1
	}
}
dev.off()
####################################################################################
mean(fstGw[,])

####################################################################################
## simple plot of Fst continuum, this includes replicates, just sorted median (across chroms) Fst
plot(sort(apply(fstGw,1,median)),pch=19)
## pretty sure the uptick at the far right is Alaskc and France
## same thing just Z
plot(sort(fstGw[,24]),pch=19) ## very similar
## and against each other
plot(apply(fstGw,1,median),fstGw[,24],pch=19,xlab="Genome median",ylab="Z chrom")
abline(a=0,b=1)
## mostly correlated, Z higher


### window examples

pdf("FstWinsABM20_SIN10.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-1;j<-28 
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()

######### plot Fst with 200 window for different pairs ########
k<-9

#### pop pair 1 vs 2: Lyc-ABM20 x Lyc-BCR17, Fst= 0.054 ####
i=1
j=2

pdf("AF_ABM20_BCR17.pdf",width=9,height=9)
par(mfrow=c(2,2))
par(mar=c(4.5,5,2.5,.5))
for(k in c(1,2,3,24)){
  
  keep<-which(n[[mord[k]]][,i] > 10 & n[[mord[k]]][,j] > 10)
  Nw<-floor(length(keep)/200)
  win<-rep(1:Nw,each=200)
  P1<-P[[mord[k]]][keep,i]
  P2<-P[[mord[k]]][keep,j]
  plot(P1,P2,pch=16,cex=0.3,xlab="Lyc-ABM20",ylab="Lyc-BCR17",
       cex.lab=1.3,main=paste("L",k),cex.main=2)
}

dev.off()
#### pop pair 1 vs 3: Lyc-ABM20 x Lyc-BHP19, Fst= 0.136 ###
i=1
j=3

pdf("AF_ABM20_BHP19.pdf",width=9,height=9)
par(mfrow=c(2,2))
par(mar=c(4.5,5,2.5,.5))
for(k in c(1,2,3,24)){
  
  keep<-which(n[[mord[k]]][,i] > 10 & n[[mord[k]]][,j] > 10)
  Nw<-floor(length(keep)/200)
  win<-rep(1:Nw,each=200)
  P1<-P[[mord[k]]][keep,i]
  P2<-P[[mord[k]]][keep,j]
  plot(P1,P2,pch=16,cex=0.3,xlab="Lyc-ABM20",ylab="Lyc-BHP19",
       cex.lab=1.3,main=paste("L",k),cex.main=2)
}

dev.off()
#### pop pair 1 vs 15: Lyc-ABM20 x Lyc-MEN12 Fst = 0.413

i=1
j=17

pdf("AF_ABM20_MEN12.pdf",width=9,height=9)
par(mfrow=c(2,2))
par(mar=c(4.5,5,2.5,.5))
for(k in c(1,2,3,24)){
  
  keep<-which(n[[mord[k]]][,i] > 10 & n[[mord[k]]][,j] > 10)
  Nw<-floor(length(keep)/200)
  win<-rep(1:Nw,each=200)
  P1<-P[[mord[k]]][keep,i]
  P2<-P[[mord[k]]][keep,j]
  plot(P1,P2,pch=16,cex=0.3,xlab="Lyc-ABM20",ylab="Lyc-MEN12",
       cex.lab=1.3,main=paste("L",k),cex.main=2)
}

dev.off()


for(k in mord){
  keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
  Nw<-floor(length(keep)/200)
  win<-rep(1:Nw,each=200)
  
  P1<-P[[k]][keep,i]
  P2<-P[[k]][keep,j]
  
  pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
  Ht<-2*pbar*(1-pbar)
  Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
  P1win<-tapply(X=P1[1:(Nw*200)],INDEX=win,mean)
  P2win<-tapply(X=P2[1:(Nw*200)],INDEX=win,mean)
  
  P1<-P[[k]][keep,i]
  P2<-P[[k]][keep,j]
  
  plot(P1,P2,pch=16,cex=0.3,xlab="Lyc-ABM20",ylab="Lyc-MEN12")
  
  plot(P1win,P2win,pch=16,cex=0.3,xlab="Lyc-ABM20",ylab="Lyc-MEN12")
}
  keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
  Nw<-floor(length(keep)/200)
  win<-rep(1:Nw,each=200)
  
  pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
  Ht<-2*pbar*(1-pbar)
  Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
  Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
  Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
  plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
  title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}


###################################################
pdf("FstWinsABM20_GNP17.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-1;j<-13
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()

pdf("FstWinsMR20_YG20.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-18;j<-27
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()

pdf("FstWinsMR20_CP19.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-18;j<-9
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}

dev.off()

pdf("FstWinsBCR17_BTB17.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-2;j<-6
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()

pdf("FstWinsMEN12_VE20.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-17;j<-26 
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()

pdf("FstWinsSBW18_VE20.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-20;j<-26 
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()
