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

#one way anova
simp.mod<-lm(LC50.PBO~nAChR.81,data=sumDat)
summary(simp.mod)
#post hoc tests
TukeyHSD(aov(simp.mod))
plot(simp.mod,1) #seems there is different variance between group
#testing heteroscedasticity
leveneTest(LC50.PBO~nAChR.81,data=sumDat) #there is heteroscedasticity

#log transformation to correct for heteroscedasticity
leveneTest(log(LC50.PBO)~nAChR.81,data=sumDat) #ok now
log.mod<-lm(log(LC50.PBO)~nAChR.81,data=sumDat)
anova(log.mod) #the group has an effect on the LC50.PBO
summary(log.mod)
#post hoc tests
TukeyHSD(aov(log.mod))
plot(log.mod,1)
#Means and SE obtained by back-transformation via the delta method
BackTrans<-emmeans(log.mod,~nAChR.81,type="response")
BackTrans


##############################################################################/
#Additional testing the different in LC50 with PBO####
##############################################################################/

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

pdf(file="output/Figure_X_TSRcomp.pdf",width=6,height=5.5)
op<-par(mar=c(5.1,5.1,1.1,1.1))
vioplot(log(sumDat$LC50.PBO)~sumDat$nAChR.81,log="",
        col=c("green3","orange3","red3"),las=1,ann=FALSE,
        bty="l",axes=FALSE,yaxt="n",xaxt="n",lwd=3,
        border=c("green4","orange4","red4"))
axis(1,at=c(1,2,3),labels=c("[RR]","[RT]","[TT]"),lwd=3,font=2)
axis(2,at=c(log(50),log(100),log(200),log(500),log(1000),log(2000)),
     labels=c("50","100","200","500","1000","2000"),lwd=3,
     las=1,font=2)
title(xlab="Target-site resistance genotype",
      ylab="LC50 with PBO",
      cex.lab=1.5,font=2)
box(bty="o",lwd=3)
stripchart(log(sumDat$LC50.PBO)~sumDat$nAChR.81,vertical=TRUE,
           method="jitter",pch=21,add=TRUE,col="grey50",
           at=c(1:3),cex=1.2,
           bg="grey80")
par(op)
dev.off()


##############################################################################/
#END
##############################################################################/