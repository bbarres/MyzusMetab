##############################################################################/
##############################################################################/
#Regression of CYP6CY3 expression against metabolic resistant proportion
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
source("MyMetab_load.R")
#limiting the dataset to the R81T-[RR] genotypes
sumDatRR<-sumDat[sumDat$nAChR.81=="[RR]",]


##############################################################################/
#Fitting different models with a maximum threshold value####
##############################################################################/

#Asymptotic regression model with 3 parameters
AssyReg.m3<-drm(sumDatRR$propMeta~sumDatRR$CY3_EXP,
                data=sumDatRR,fct=AR.3())
plot(AssyReg.m3,type="confidence",log="",col="blue",lwd=3,lty=2)
plot(AssyReg.m3,type="obs",add=TRUE)
summary(AssyReg.m3)

#Asymptotic regression model with 2 parameters
AssyReg.m2<-drm(sumDatRR$propMeta~sumDatRR$CY3_EXP,
                data=sumDatRR,fct=AR.2())
plot(AssyReg.m2,type="confidence",log="",col="blue",lwd=3,lty=2)
plot(AssyReg.m2,type="obs",add=TRUE)
summary(AssyReg.m2)

anova(AssyReg.m2,AssyReg.m3) #the model with fewer parameters is preferable


#Michaelis-Menten model with 3 parameters
Micha.m3<-drm(sumDatRR$propMeta~sumDatRR$CY3_EXP,
              data=sumDatRR,fct=MM.3())
plot(Micha.m3,type="confidence",log="",col="blue",lwd=3,lty=2)
plot(Micha.m3,type="obs",add=TRUE)
summary(Micha.m3)

#Michaelis-Menten model with 2 parameters
Micha.m2<-drm(sumDatRR$propMeta~sumDatRR$CY3_EXP,
              data=sumDatRR,fct=MM.2())
plot(Micha.m2,type="confidence",log="",col="blue",lwd=3,lty=2)
plot(Micha.m2,type="obs",add=TRUE)
summary(Micha.m2)

anova(Micha.m2,Micha.m3) #the model with fewer parameters is preferable

#results of the final model, the one with the smallest RSE: Micha.m2
summary(Micha.m2)
coef(Micha.m2)


##############################################################################/
#Figure 4: Plotting the asymptotic regression model####
##############################################################################/

modplot<-Micha.m2
colovec=c("red3","orange3","green3")

pdf(file="output/Figure_4_metaRR.pdf",width=7,height=6)
op<-par(mar=c(5.1,5.1,1.1,1.1))
plot(modplot,type="confidence",log="",col="skyblue3",
     lwd=4,lty=2,ann=FALSE,axes=FALSE,bty="n",
     xlab="",
     ylab="",
     ylim=c(0,1),xlim=c(0,24))
points(sumDat$propMeta~sumDat$CY3_EXP,
       bg=colovec[as.numeric(sumDat$nAChR.81)],pch=21)
plot(modplot,type="obs",pch=19,col="green3",cex=1.5,add=TRUE)
axis(1,lwd=3,font=2)
axis(2,lwd=3,las=1,font=2)
box(bty="o",lwd=3)
title(xlab="CYP6CY3 expression level",
      ylab="Proportion of metabolic resistance",
      cex.lab=1.4,font=2,line=3.5)
legend(18,0.5,pch=c(21,21,19),pt.cex=c(1.5,1.5,2),
       pt.bg=colovec,col=c("black","black","green3"),
       legend=c(expression("81"^"TT"),
                expression("81"^"RT"),
                expression("81"^"RR")),
       y.intersp=1.4,bty="n")
par(op)
dev.off()


##############################################################################/
#END
##############################################################################/

##############################################################################/
#Additional code for Negative exponential model fitting####
##############################################################################/

#adapted from https://rpubs.com/mengxu/exponential-model
#Prepare a good inital state
theta.0 <- max(sumDatRR$propMeta) * 1.1
model.0 <- lm(log(-propMeta+theta.0)~CY3_EXP,data=sumDatRR)
alpha.0 <- -exp(coef(model.0)[1])
beta.0 <- coef(model.0)[2]

debut<-list(alpha=alpha.0,beta=beta.0,theta=theta.0)

#Fit the model
negExpo.m1<-nls(propMeta~alpha*exp(beta*CY3_EXP)+theta,
                data=sumDatRR,start=debut)

#add fitted curve
plot(sumDatRR$CY3_EXP,sumDatRR$propMeta)
lines(seq(0.01,20.5,by=0.05),
      predict(negExpo.m1,list(CY3_EXP=seq(0.01,20.5,by=0.05)),
              interval="confidence",level=0.95),
      col='darkgreen',lwd=3)
summary(negExpo.m1)