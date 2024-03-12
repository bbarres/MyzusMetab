##############################################################################/
##############################################################################/
#Comparison between LD50 with or without PBO 
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
source("MyMetab_load.R")


##############################################################################/
#Correlation between LD50 with or without PBO####
##############################################################################/

#comparing the variance between LC50 and LC50.PBO
mean(sumDat$LC50)
tapply(sumDat$LC50,sumDat$nAChR.81,mean)
var(sumDat$LC50)
tapply(sumDat$LC50,sumDat$nAChR.81,var)
mean(sumDat$LC50.PBO)
tapply(sumDat$LC50.PBO,sumDat$nAChR.81,mean)
var(sumDat$LC50.PBO)
tapply(sumDat$LC50.PBO,sumDat$nAChR.81,var)

#it seems there is an important reduction in variance
var(sumDat$LC50)/var(sumDat$LC50.PBO) #almost 30 times greater for LC50 
                                      #than for LC50.PBO
tapply(sumDat$LC50,sumDat$nAChR.81,var)/
  tapply(sumDat$LC50.PBO,sumDat$nAChR.81,var)
#testing the difference between variance using a Levene's test
long_sumDat<-pivot_longer(sumDat[,1:5],cols=4:5,
                   names_to="treatment",values_to="LC50value")
leveneTest(LC50value~treatment,data=long_sumDat,center="mean")

#investigating the correlation between thiacloprid LC50 with or without PBO
withwith.mod<-lm(sumDat$LC50.PBO~sumDat$LC50)
#the residuals are not normally distributed
plot(withwith.mod,1)
plot(withwith.mod,2)
#not surprising since we have seen the variance significantly differs
leveneTest(LC50value~treatment,data=long_sumDat,center="mean")
#we log transform both the LC50 values to fix the variance issue
leveneTest(log(LC50value)~treatment,data=long_sumDat,center="mean")
withwith.mod<-lm(log(sumDat$LC50.PBO)~log(sumDat$LC50))
#the residuals look better
plot(withwith.mod,1)
plot(withwith.mod,2)

#the results of the correlation between LC50 and LC50.PBO values
summary(withwith.mod)
anova(withwith.mod)


##############################################################################/
#Figure 1A: Correlation between LC50 with or without PBO####
##############################################################################/

#producing the plot
pdf(file="output/Figure_1A_CorrWithWitho.pdf",width=6,height=5.5)
op<-par(mar=c(5.1,5.1,1.1,1.1))
colovec=c("red3","orange3","green3")
plot(log(sumDat$LC50.PBO)~log(sumDat$LC50),log="",pch=21,
     bg=colovec[as.numeric(sumDat$nAChR.81)],las=1,
     ann=FALSE,axes=FALSE,cex=2)
abline(lm(log(sumDat$LC50.PBO)~log(sumDat$LC50)),
       col="black",lwd=4,lty=2,untf=TRUE)
points(log(sumDat$LC50.PBO)~log(sumDat$LC50),pch=21,cex=2,
       bg=colovec[as.numeric(sumDat$nAChR.81)])
axis(1,at=log(c(50,100,200,500,1000,2500,7500,20000)),lwd=3,
     labels=c("50","100","200","500","1000","2500","7500","20000"),
     font=2)
axis(2,at=log(c(50,100,200,500,1000,2000)),
     labels=c("50","100","200","500","1000","2000"),lwd=3,
     las=1,font=2)
legend(log(50),log(2500),pch=21,pt.cex=2,pt.bg=colovec,
       legend=c(expression("81"^"TT"),
                expression("81"^"RT"),
                expression("81"^"RR")),
       y.intersp=1.2,
       bty="n")
title(xlab=expression(paste("Thiacloprid LC50 without PBO (",mu,"g/L)")),
      ylab=expression(paste("Thiacloprid LC50 with PBO (",mu,"g/L)")),
      cex.lab=1.4,font=2,line=3.5)
box(bty="o",lwd=3)
par(op)
dev.off()


##############################################################################/
#Figure 1B: LC50 dumbellplot with or without PBO####
##############################################################################/

#this code was adapted from: 
#https://r-coder.com/dot-plot-r/?utm_content=cmp-true

pdf(file="output/Figure_1B_LC50dumbell.pdf",width=9,height=7)
op<-par(mar=c(5.1,4.1,1.1,1.1))
colovec=c("red3","orange3","green3")
#ordering the data by increasing LC50
sumDat<-sumDat[order(sumDat$nAChR.81,-sumDat$LC50,decreasing=TRUE),]

#ordPlot<-sort.list(as.numeric(sumDat$nAChR.81),decreasing=TRUE)
ycorec<-cumsum(c(0,diff(as.numeric(sumDat$nAChR.81)) != 0))
ycoor<-1:dim(sumDat)[1] + 2*ycorec

dotchart(sumDat$LC50,
         groups=sumDat$nAChR.81,
         labels=sumDat$clone.ID,
         xlim=range(sumDat$LC50.PBO,sumDat$LC50)+c(-10,1000),
         col=colovec[as.numeric(sumDat$nAChR.81)],
         log="x",pch="")

for(i in 1:dim(sumDat)[1]) {
  segments(min(sumDat$LC50.PBO[i],sumDat$LC50[i]),ycoor[i],
           max(sumDat$LC50.PBO[i],sumDat$LC50[i]),ycoor[i],
           lwd=5,col=colovec[as.numeric(sumDat$nAChR.81)][i]) 
}

points(sumDat$LC50,ycoor,pch=19,cex=2)
points(sumDat$LC50.PBO,ycoor,pch=21,cex=2,bg="white")

legend(2500,10,legend=c("without PBO","with PBO"),bg="white",
       pch=c(19,21),col=c("black","black"),pt.bg=c("black","white"),
       title="Bioassays:",title.font=2,title.cex=1.2,pt.cex=2,
       y.intersp=1.2)

title(xlab=expression(paste("Thiacloprid LC50 (",mu,"g/L)")),
      ylab="",
      cex.lab=1.5,font=2)

par(op)
dev.off()


##############################################################################/
#END
##############################################################################/