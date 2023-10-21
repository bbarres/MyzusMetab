##############################################################################/
##############################################################################/
#Impact of target-site genotype on the resistance phenotype 
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
source("MyMetab_load.R")


##############################################################################/
#Testing the different in LC50 with PBO####
##############################################################################/

vioplot(log(sumDat$LC50.PBO)~sumDat$nAChR.81,log="",
        col=c("green3","orange3","red3"),las=1)
stripchart(log(sumDat$LC50.PBO)~sumDat$nAChR.81,vertical=TRUE,
           method="jitter",pch=21,add=TRUE,col="black",
           at=c(1:3),cex=0.8)

#one way anova
simp.mod<-lm(LC50.PBO~nAChR.81,data=sumDat)
summary(simp.mod)
#post hoc tests
TukeyHSD(aov(simp.mod))
plot(simp.mod,1) #seems there is different variance between group
#testing heteroscedasticity
leveneTest(LC50.PBO~nAChR.81,data=sumDat) #there is heteroscedasticity

#log transformation to correct for heteroscedasticity
log.mod<-lm(log(LC50.PBO)~nAChR.81,data=sumDat)
summary(log.mod)
#post hoc tests
TukeyHSD(aov(log.mod))
plot(log.mod,1)
#Means and SE obtained by back-transformation via the delta method
BackTrans<-emmeans(log.mod,~nAChR.81,type="response")


#anova with Welch correction
oneway.test(sumDat$LC50.PBO~sumDat$nAChR.81,var.equal=FALSE)
pairwise.t.test(sumDat$LC50.PBO,sumDat$nAChR.81,
                pool.sd=FALSE,p.adjust.method="BH")

#another more conservative methods: non parametric test
kruskal.test((sumDat$LC50.PBO)~sumDat$nAChR.81)
#pairwise test
pairwise.wilcox.test(sumDat$LC50.PBO,sumDat$nAChR.81,p.adjust.method="BH")


##############################################################################/
#Figure XX: plot of the distribution of LC50 with PBO by TSR genotype####
##############################################################################/

vioplot(log(sumDat$LC50.PBO)~sumDat$nAChR.81,log="",
        col=c("green3","orange3","red3"),las=1,ann=FALSE,
        bty="l",axes=FALSE,yaxt="n")
box(bty="l")

stripchart(log(sumDat$LC50.PBO)~sumDat$nAChR.81,vertical=TRUE,
           method="jitter",pch=21,add=TRUE,col="black",
           at=c(1:3),cex=0.8)


##############################################################################/
#END
##############################################################################/