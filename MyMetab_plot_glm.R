##############################################################################/
##############################################################################/
#plot probablement à jeter
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
source("MyMetab_glm.R")


##############################################################################/
#Regression analyses####
##############################################################################/

#figure for the mod1
plot(sumDat$CY3_EXP,log(sumDat$LC50),type="n",las=1)
points(sumDat[sumDat$nAChR.81=="[TT]",]$CY3_EXP,
       log(sumDat[sumDat$nAChR.81=="[TT]",]$LC50),pch=16,col="red")
points(sumDat[sumDat$nAChR.81=="[RT]",]$CY3_EXP,
       log(sumDat[sumDat$nAChR.81=="[RT]",]$LC50),pch=16,col="orange")
points(sumDat[sumDat$nAChR.81=="[RR]",]$CY3_EXP,
       log(sumDat[sumDat$nAChR.81=="[RR]",]$LC50),pch=16,col="green")

genoAch<-factor(rep("[TT]",251))
xv<-seq(0,25,0.1)
yv<-predict(mod1,list(CY3_EXP=xv,nAChR.81=genoAch))
lines(xv,yv,col="red")

genoAch<-factor(rep("[RT]",251))
xv<-seq(0,25,0.1)
yv<-predict(mod1,list(CY3_EXP=xv,nAChR.81=genoAch))
lines(xv,yv,col="orange")

genoAch<-factor(rep("[RR]",251))
xv<-seq(0,25,0.1)
yv<-predict(mod1,list(CY3_EXP=xv,nAChR.81=genoAch))
lines(xv,yv,col="green")


##############################################################################/
#END
##############################################################################/