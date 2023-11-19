##############################################################################/
##############################################################################/
#Modelling the LD50
##############################################################################/
##############################################################################/


#loading the packages necessary for the analysis
source("MyMetab_load.R")
#we reorder the level of the nAChR.81 factor so that the sensitive genotype
#is the reference
levels(sumDat$nAChR.81)<-c("[RR]","[RT]","[TT]")

##############################################################################/
#Correlation between P450 quantitative variables####
##############################################################################/

str(sumDat)

#a function to compute the absolute correlation between pairs of variables
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use="pairwise.complete.obs"))
  txt <- format(c(r, 2), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

#because we suspect collinearity between genes copy number, we 
#first look at correlation between the different potential variable
pairs(sumDat[,c(6:11)],upper.panel=panel.cor,las=1)

#producing the supplementary figure
pdf(file="output/Figure_SX_pairs.pdf",width=9,height=7)
colovec=c("red3","orange3","green3")
pairs(sumDat[,c(6:11)],upper.panel=panel.cor,las=1,
      col=colovec[as.numeric(sumDat$nAChR.81)],pch=19)
dev.off()

#because the correlation between copy number of the different genes
#and level of expression of the different genes are correlated, we 
#pick only one of these variables. Since CY23 is a less good candidate
#we don't include the related variables either


##############################################################################/
#Actual modeling####
##############################################################################/

#the complete model including the R81T genotype, one P450 quantification 
#variable, the genetic group and their interactions. The different P450 
#quantification variables are tested to select the best one
modT<-glm(LC50~nAChR.81*CY3_EXP*genetic.group,data=sumDat,
          family=stats::gaussian(link="log"))
summary(modT)$aic #AIC 324.388
modT<-glm(LC50~nAChR.81*CY3_CN*genetic.group,data=sumDat,
          family=stats::gaussian(link="log"))
summary(modT)$aic #AIC 427.2367
modT<-glm(LC50~nAChR.81*CY4_EXP*genetic.group,data=sumDat,
          family=stats::gaussian(link="log"))
summary(modT)$aic #AIC 394.7389
modT<-glm(LC50~nAChR.81*CY4_CN*genetic.group,data=sumDat,
          family=stats::gaussian(link="log"))
summary(modT)$aic #AIC 417.6384
#the model with the smallest AIC is the one using CY3_EXP variable
#therefore we select this variable for further analyses


modT<-glm(LC50~genetic.group*nAChR.81*CY3_EXP,data=sumDat,
          family=stats::gaussian(link="log"))
summary(modT)
#the 3 factor interaction is irrelevant here, so we remove it
modT<-glm(LC50~nAChR.81+CY3_EXP+genetic.group+
            genetic.group:nAChR.81+genetic.group:CY3_EXP+nAChR.81:CY3_EXP,
          data=sumDat,family=stats::gaussian(link="log"))
summary(modT)
anova(modT,test="Chisq")
#the interaction between genetic.group and other factors are not significant
#so we remove them too
modT<-glm(LC50~nAChR.81+CY3_EXP+nAChR.81:CY3_EXP+genetic.group,
          data=sumDat,family=stats::gaussian(link="log"))
anova(modT,test="Chisq")
summary(modT)
#the genetic.group doesn't seem to have an impact on LC50, so we remove this 
#variable from the model
mod1<-glm(LC50~nAChR.81*CY3_EXP,data=sumDat,
          family=stats::gaussian(link="log"))
summary(mod1)$aic #the AIC is improved
#let's try to remove the interaction to simplify further the model
mod2<-glm(LC50~nAChR.81+CY3_EXP,data=sumDat,
          family=stats::gaussian(link="log"))
summary(mod2)$aic #the AIC is greater than the one obtained with mod1
anova(mod2,mod1,test="Chisq") #the model with the interaction is better

#here are the results for the final model with the nAChR.81 genotype, 
#the CY3 expression level and their interaction
summary(mod1)
anova(mod1,test="Chisq")
plot(mod1,1)


##############################################################################/
#END
##############################################################################/


#modelling the EC50 with PBO
modPBO.1<-glm(LC50.PBO~nAChR.81*CY3_EXP,data=sumDat,
              family=stats::gaussian(link="log"))
summary(modPBO.1) #best model AIC 282
anova(modPBO.1,test="Chisq")
plot(modPBO.1)
modPBO.2<-glm(LC50.PBO~nAChR.81*CY4_EXP,data=sumDat,
              family=stats::gaussian(link="log"))
summary(modPBO.2) #AIC 324
modPBO.3<-glm(LC50.PBO~nAChR.81+CY3_EXP,data=sumDat,
              family=stats::gaussian(link="log"))
summary(modPBO.3) #AIC 292
modPBO.4<-glm(LC50.PBO~nAChR.81+CY4_EXP,data=sumDat,
              family=stats::gaussian(link="log"))
summary(modPBO.4) #AIC 322
anova(modPBO.3,modPBO.1,test="Chisq")

#modelling the difference/ratio with or without PBO
modDif1<-glm(LC50/LC50.PBO~nAChR.81*CY3_EXP,data=sumDat,
             family=stats::gaussian(link="identity"))
summary(modDif1) #best model AIC 76
anova(modDif1,test="Chisq")
modDif2<-glm(LC50/LC50.PBO~nAChR.81*CY4_EXP,data=sumDat,
             family=stats::gaussian(link="identity"))
summary(modDif2) #AIC 77
modDif3<-glm(LC50/LC50.PBO~nAChR.81+CY3_EXP,data=sumDat,
             family=stats::gaussian(link="identity"))
summary(modDif3) #AIC 80 
modDif4<-glm(LC50/LC50.PBO~nAChR.81+CY4_EXP,data=sumDat,
             family=stats::gaussian(link="identity"))
summary(modDif4) #AIC 81
anova(modDif3,modDif1,test="Chisq")

#to sum up the best model
summary(mod1)
summary(modPBO.1)
summary(modDif1)
