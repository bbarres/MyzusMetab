##############################################################################/
##############################################################################/
#Plot of the copy number and expression level of CYP6CY3
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
source("MyMetab_load.R")
#we reorder the level of the nAChR.81 factor so that the sensitive genotype
#is the reference
sumDat$nAChR.81<-factor(sumDat$nAChR.81,levels=c("[RR]","[RT]","[TT]"))


##############################################################################/
##############################################################################/

error.bar<-function(x,y,upper,lower=upper,length=0.03,COL,...){                         
  if(length(x)!=length(y)|
     length(y)!=length(lower)|
     length(lower)!=length(upper))   
    stop("vectors must be same length")                                                       
  arrows(x,y+upper,x,y-lower,angle=90,code=3,length=length,lwd=1,col=COL, ...)                  
}


##### Figure 2

pdf(file="Figure 2 Copy number and expression level CY3 CY4.pdf",width=7,height=5)
par(mfrow=c(2,2),mar=(c(6, 4, 2, 2) + 0.1))
espace_entre_barre<-c(0,rep(0.3,11),0.9,rep(0.3,6),0.9,rep(0.3,3))

colovec=c("green3" ,"orange3","red3")

#sort by genotype R81T and LD50
sumDat<-sumDat[do.call(order, sumDat[c(3,4)]),]

#CYP6CY3
barplot(names.arg=sumDat$clone.ID,height=sumDat$CY3_CN,xlab="",ylim=c(0,10),main="",las=2,border="black",cex.lab = 0.8,
        col=colovec[sumDat$nAChR.81],space=espace_entre_barre,axisnames=F)
coord<-barplot(plot=F,names.arg=sumDat$clone.ID,height=sumDat$CY3_CN,space=espace_entre_barre)
error.bar(coord,sumDat$CY3_CN, sumDat$CY3_SE_CN,COL="black")
axis(1,at=c(-5,coord),labels=rep("",length(coord)+1))
title(ylab="CY3 relative copy number", line=2.5, cex.lab=0.8)
mtext(side=1,at=c(-5,coord),line=4,text=c("",as.character(sumDat$clone.ID)),col=c("white",colovec[sumDat$nAChR.81]),las=2,cex=0.5,adj=0)
mtext(side=1,at=c(8,21,28.5),line=4.6,text=expression(81^RR,81^RT,81^TT),cex=0.9)
mtext(side=3,at=1.5,line=0,text="A",cex=1.7)

#CY4
barplot(names.arg=sumDat$clone.ID,height=sumDat$CY4_CN,xlab="",ylim=c(0,10),main="",las=2,border="black",cex.lab = 0.8,
        col=colovec[sumDat$nAChR.81],cex.names=0.7,space=espace_entre_barre, axisnames=F)
coord<-barplot(plot=F,names.arg=sumDat$clone.ID,height=sumDat$CY4_CN,space=espace_entre_barre)
error.bar(coord,sumDat$CY4_CN,sumDat$CY4_SE_CN,COL="black")
axis(1,at=c(-5,coord),labels=rep("",length(coord)+1))
title(ylab="CY4 relative copy number", line=2.5, cex.lab=0.8)
mtext(side=1,at=c(-5,coord),line=4,text=c("",as.character(sumDat$clone.ID)),col=c("white",colovec[sumDat$nAChR.81]),las=2,cex=0.5,adj=0)
mtext(side=1,at=c(8,21,28.5),line=4.6,text=expression(81^RR,81^RT,81^TT),cex=0.9)
mtext(side=3,at=1.5,line=0,text="B",cex=1.7)

#CY3 ARN
barplot(names.arg=sumDat$clone.ID,height=sumDat$CY3_EXP,xlab="",ylim=c(0,33),main="",las=2,border="black",cex.lab = 0.8,
        col=colovec[sumDat$nAChR.81],space=espace_entre_barre, axisnames=F)
coord<-barplot(plot=F,names.arg=sumDat$clone.ID,height=sumDat$CY3_EXP,space=espace_entre_barre)
error.bar(coord,sumDat$CY3_EXP, sumDat$CY3_SE_EXP,COL="black")
axis(1,at=c(-5,coord),labels=rep("",length(coord)+1))
title(ylab="CY3 relative expression level", line=2.5, cex.lab=0.8)
mtext(side=1,at=c(-5,coord),line=4,text=c("",as.character(sumDat$clone.ID)),col=c("white",colovec[sumDat$nAChR.81]),las=2,cex=0.5,adj=0)
mtext(side=1,at=c(8,21,28.5),line=4.6,text=expression(81^RR,81^RT,81^TT),cex=0.9)
mtext(side=3,at=1.5,line=-0.4,text="C",cex=1.7)

#CY4 ARN
barplot(names.arg=sumDat$clone.ID,height=sumDat$CY4_EXP,xlab="",ylim=c(0,33),main="",las=2,border="black",cex.lab = 0.8,
        col=colovec[sumDat$nAChR.81],cex.names=0.7,space=espace_entre_barre, axisnames=F)
coord<-barplot(plot=F,names.arg=sumDat$clone.ID,height=sumDat$CY4_EXP,space=espace_entre_barre)
error.bar(coord,sumDat$CY4_EXP,sumDat$CY4_SE_EXP,COL="black")
axis(1,at=c(-5,coord),labels=rep("",length(coord)+1))
title(ylab="CY4 relative expression level", line=2.5, cex.lab=0.8)
mtext(side=1,at=c(-5,coord),line=4,text=c("",as.character(sumDat$clone.ID)),col=c("white",colovec[sumDat$nAChR.81]),las=2,cex=0.5,adj=0)
mtext(side=1,at=c(8,21,28.5),line=4.6,text=expression(81^RR,81^RT,81^TT),cex=0.9)
mtext(side=3,at=1.5,line=-0.4,text="D",cex=1.7)

dev.off()




##### Figure S1

pdf(file="Figure S1 Copy number and expression level CY23.pdf",width=4,height=6)
par(mfrow=c(2,1),mar=(c(6, 4, 1, 2) + 0.1))
espace_entre_barre<-c(0,rep(0.3,11),0.9,rep(0.3,6),0.9,rep(0.3,3))

colovec=c("green3" ,"orange3","red3")

#CY23
barplot(names.arg=sumDat$clone.ID,height=sumDat$CY23_CN,xlab="",ylim=c(0,10),ylab="CY23 relative copy number",main="",las=2,border="black",cex.lab = 0.8,
        col=colovec[sumDat$nAChR.81],cex.names=0.7,space=espace_entre_barre, axisnames=F)
coord<-barplot(plot=F,names.arg=sumDat$clone.ID,height=sumDat$CY23_CN,space=espace_entre_barre)
error.bar(coord,sumDat$CY23_CN,sumDat$CY23_SE_CN,COL="black")
axis(1,at=c(-5,coord),labels=rep("",length(coord)+1))
mtext(side=1,at=c(-5,coord),line=3.6,text=c("",as.character(sumDat$clone.ID)),col=c("white",colovec[sumDat$nAChR.81]),las=2,cex=0.5,adj=0)
mtext(side=1,at=c(8,21,28.5),line=4.2,text=expression(81^RR,81^RT,81^TT),cex=0.9)
mtext(side=3,at=1.5,line=-0.4,text="A",cex=1.7)

#CY23 ARN
barplot(names.arg=sumDat$clone.ID,height=sumDat$CY23_EXP,xlab="",ylim=c(0,33),ylab="CY23 relative expression level",main="",las=2,border="black",cex.lab = 0.8,
        col=colovec[sumDat$nAChR.81],cex.names=0.7,space=espace_entre_barre, axisnames=F)
coord<-barplot(plot=F,names.arg=sumDat$clone.ID,height=sumDat$CY23_EXP,space=espace_entre_barre)
error.bar(coord,sumDat$CY23_EXP,sumDat$CY23_SE_EXP,COL="black")
axis(1,at=c(-5,coord),labels=rep("",length(coord)+1))
mtext(side=1,at=c(-5,coord),line=3.6,text=c("",as.character(sumDat$clone.ID)),col=c("white",colovec[sumDat$nAChR.81]),las=2,cex=0.5,adj=0)
mtext(side=1,at=c(8,21,28.5),line=4.2,text=expression(81^RR,81^RT,81^TT),cex=0.9)
mtext(side=3,at=1.5,line=-0.4,text="B",cex=1.7)

dev.off()




