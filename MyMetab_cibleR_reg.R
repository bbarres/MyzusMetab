##############################################################################/
##############################################################################/
#Impact of target-site genotype on the resistance phenotype 
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
source("MyMetab_load.R")
#we reorder the level of the nAChR.81 factor so that the sensitive genotype
#is the reference
sumDat$nAChR.81<-factor(sumDat$nAChR.81,levels=c("[RR]","[RT]","[TT]"))


##############################################################################/
#Testing the different in LC50 between R81T genotype groups####
##############################################################################/

#one way anova
simp.mod<-lm(LC50~nAChR.81,data=sumDat)
summary(simp.mod)
#post hoc tests
TukeyHSD(aov(simp.mod))
plot(simp.mod,1) #seems there is different variance between group
#testing heteroscedasticity
leveneTest(LC50~nAChR.81,data=sumDat,
           center="mean") #there is heteroscedasticity

#log transformation to correct for heteroscedasticity
leveneTest(log(LC50)~nAChR.81,data=sumDat,center="mean") #ok now
log.mod<-lm(log(LC50)~nAChR.81,data=sumDat)
anova(log.mod) #the group has an effect on the LC50.PBO
summary(log.mod)
#post hoc tests
TukeyHSD(aov(log.mod))
plot(log.mod,1)
plot(log.mod,2)
exp(TukeyHSD(aov(log.mod))$nAChR.81[,1])
exp(TukeyHSD(aov(log.mod))$nAChR.81[,2])
exp(TukeyHSD(aov(log.mod))$nAChR.81[,3])

#Means and SE obtained by back-transformation via the delta method
BackTrans<-emmeans(log.mod,~nAChR.81,type="response")
BackTrans


##############################################################################/
#Testing the different in LC50 with PBO####
##############################################################################/

#one way anova
simp.modPBO<-lm(LC50.PBO~nAChR.81,data=sumDat)
summary(simp.modPBO)
#post hoc tests
TukeyHSD(aov(simp.modPBO))
plot(simp.modPBO,1) #seems there is different variance between group
#testing heteroscedasticity
leveneTest(LC50.PBO~nAChR.81,data=sumDat,
           center="mean") #there is heteroscedasticity

#log transformation to correct for heteroscedasticity
leveneTest(log(LC50.PBO)~nAChR.81,data=sumDat,center="mean") #ok now
log.modPBO<-lm(log(LC50.PBO)~nAChR.81,data=sumDat)
anova(log.modPBO) #the group has an effect on the LC50.PBO
summary(log.modPBO)
#post hoc tests
TukeyHSD(aov(log.modPBO))
plot(log.modPBO,1)
plot(log.modPBO,2)
exp(TukeyHSD(aov(log.modPBO))$nAChR.81[,1])
exp(TukeyHSD(aov(log.modPBO))$nAChR.81[,2])
exp(TukeyHSD(aov(log.modPBO))$nAChR.81[,3])

#Means and SE obtained by back-transformation via the delta method
BackTransPBO<-emmeans(log.modPBO,~nAChR.81,type="response")
BackTransPBO


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
#Figure 3: plot of the distribution of LC50 with PBO by TSR genotype####
##############################################################################/

pdf(file="output/Figure_3_TSRcomp.pdf",width=6,height=5.5)
op<-par(mar=c(5.1,5.1,1.1,1.1))
vioplot(log(sumDat$LC50.PBO)~sumDat$nAChR.81,log="",
        col=c("green3","orange3","red3"),las=1,ann=FALSE,
        bty="l",axes=FALSE,yaxt="n",xaxt="n",lwd=3,wex=0.8,
        border=c("green4","orange4","red4"))
axis(1,at=c(1,2,3),labels=c(expression("81"^"RR"),
                            expression("81"^"RT"),
                            expression("81"^"TT")),
     lwd=3,font=2)
axis(2,at=c(log(50),log(100),log(200),log(500),log(1000),log(2000)),
     labels=c("50","100","200","500","1000","2000"),lwd=3,
     las=1,font=2)
title(xlab="Target-site resistance genotype",
      ylab=expression(paste("LC50 with PBO (",mu,"g/L)")),
      cex.lab=1.4,font=2,line=3.5)
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