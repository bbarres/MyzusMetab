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
#END
##############################################################################/