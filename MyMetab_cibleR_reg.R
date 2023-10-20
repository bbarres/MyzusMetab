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

vioplot(log(sumDat$LC50.PBO)~sumDat$nAChR.81,log="",
        col=c("green3","orange3","red3"),las=1)
stripchart(log(sumDat$LC50.PBO)~sumDat$nAChR.81,vertical=TRUE,
           method="jitter",pch=21,add=TRUE,col="black",
           at=c(1:3),cex=0.8)
TukeyHSD(aov(log(sumDat$LC50.PBO)~sumDat$nAChR.81))
temp<-(lm(log(sumDat$LC50)~sumDat$nAChR.81))


##############################################################################/
#END
##############################################################################/