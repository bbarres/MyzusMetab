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
modplot<-Micha.m2


##############################################################################/
#Assessing the equivalent model with CYP6CY4 expression covariate####
##############################################################################/

#Michaelis-Menten model with 2 parameters
Micha.m2<-drm(sumDatRR$propMeta~sumDatRR$CY4_EXP,
              data=sumDatRR,fct=MM.2())
plot(Micha.m2,type="confidence",log="",col="blue",lwd=3,lty=2)
plot(Micha.m2,type="obs",add=TRUE)
summary(Micha.m2)


##############################################################################/
#END
##############################################################################/