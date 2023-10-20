##############################################################################/
##############################################################################/
#Regression of CYP6CY3 expression against metabolic resistant proportion
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
source("MyMetab_load.R")
#limiting the dataset to the R81T-[RR] genotypes
sumDatRR<-sumDat[sumDat$nAChR.81=="RR",]
#adding a column for the proportion of LD50 linked to metabolic resistance
sumDatRR$propMeta<-(1-(sumDatRR$LC50.PBO/sumDatRR$LC50))
sumDatRR<-sumDatRR[order(sumDatRR$CY3_EXP),]



##############################################################################/
#Fitting a Weibull model####
##############################################################################/

Weib.m1<-drm(sumDatRR$propMeta~sumDatRR$CY3_EXP,
             data=sumDatRR,fct=W1.2())
plot(Weib.m1,type="confidence",log="",col="blue",lwd=3,lty=2)
plot(Weib.m1,type="obs",add=TRUE)
summary(Weib.m1)

Micha.m1<-drm(sumDatRR$propMeta~sumDatRR$CY3_EXP,
             data=sumDatRR,fct=MM.2())
plot(Micha.m1)
summary(Micha.m1)

ExpoDec.m1<-drm(sumDatRR$propMeta~sumDatRR$CY3_EXP,
                data=sumDatRR,fct=EXD.3())
plot(ExpoDec.m1)
summary(ExpoDec.m1)

#Asymptotic regression model with 2 parameters
pdf(file="output/Figure_X_metaRR.pdf",width=7,height=6)
AssyReg.m1<-drm(sumDatRR$propMeta~sumDatRR$CY3_EXP,
              data=sumDatRR,fct=AR.2())
plot(AssyReg.m1,type="confidence",log="",col="skyblue3",lwd=3,lty=2,
     xlab="CYP6CY3 expression level",
     ylab="Proportion of metabolic resistance",
     ylim=c(0,1))
plot(AssyReg.m1,type="obs",pch=19,col="green3",add=TRUE)
dev.off()
summary(AssyReg.m1)


##############################################################################/
#END
##############################################################################/


##############################################################################/
#Additional code for Negative exponential model fitting####
##############################################################################/

#adapted from https://rpubs.com/mengxu/exponential-model
# Prepare a good inital state
theta.0 <- max(sumDatRR$propMeta) * 1.1
model.0 <- lm(log(-propMeta+theta.0)~CY3_EXP,data=sumDatRR)
alpha.0 <- -exp(coef(model.0)[1])
beta.0 <- coef(model.0)[2]

debut<-list(alpha=alpha.0,beta=beta.0,theta=theta.0)

# Fit the model
negExpo.m1<-nls(propMeta~alpha*exp(beta*CY3_EXP)+theta,
                data=sumDatRR,start=debut)

# add fitted curve
plot(sumDatRR$CY3_EXP,sumDatRR$propMeta)
lines(seq(0.01,20.5,by=0.05),
      predict(negExpo.m1,list(CY3_EXP=seq(0.01,20.5,by=0.05)),
              interval="confidence",level=0.95),
      col='darkgreen',lwd=3)
summary(negExpo.m1)
