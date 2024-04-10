##############################################################################/
##############################################################################/
#Modelling the LD50
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
source("MyMetab_load.R")
#we reorder the level of the nAChR.81 factor so that the sensitive genotype
#is the reference
sumDat$nAChR.81<-factor(sumDat$nAChR.81,levels=c("[RR]","[RT]","[TT]"))


##############################################################################/
#Correlation between P450 quantitative variables####
##############################################################################/

str(sumDat)

#a function to compute the absolute correlation between pairs of variables
panel.cor <- function(x,y,digits=2,prefix="",cex.cor,...)
{
  usr<-par("usr"); on.exit(par(usr))
  par(usr=c(0,1,0,1))
  r<-abs(cor(x,y,use="pairwise.complete.obs"))
  txt<-format(c(r,2),digits=digits)[1]
  txt<-paste(prefix,txt,sep="")
  if(missing(cex.cor)) cex.cor<-0.8/strwidth(txt)
  text(0.5,0.5,txt,cex=cex.cor * r)
}

#because we suspect collinearity between genes copy number, we 
#first look at correlation between the different potential variable
pairs(sumDat[,c(6:11)],upper.panel=panel.cor,las=1)

#producing the supplementary figure S3
pdf(file="output/Figure_S2_pairs.pdf",width=9,height=7)
colovec=c("green3","orange3","red3")
pairs(sumDat[,c(6:11)],upper.panel=panel.cor,las=1,col=grey(0.0,1.0),
      bg=colovec[as.numeric(sumDat$nAChR.81)],pch=21,
      cex=1.5,cex.cor=3)
dev.off()

#because the correlation between copy number of the different genes
#and level of expression of the different genes are correlated, we 
#need to pick only one of these variables to include in the model. 


##############################################################################/
#Choice of the P450 variable####
##############################################################################/

#In order to select the most relevant P450 variable, we compare the complete
#model with the different possible variable. Since CY23 was not involved in 
#the resistance to neonicotinoid, the two related variable (number of copy 
#and expression level were not tested). The "best" P450 variable was selected
#using AIC of the full models (the smaller AIC the better)
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
#therefore we select this variable for further modelling


##############################################################################/
#Model selection using a backward stepwise regression approach####
##############################################################################/

#First we fit the complete model (including all the interaction)
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
#plot to compare the two models
plot_summs(mod1,mod2,plot.distributions=TRUE,
           model.names=c("with interactions","without interactions"))
plot_summs(mod1,mod2,plot.distributions=FALSE,
           model.names=c("with interactions","without interactions"))
#another more straightforward way to do backward elimination
step(modT,direction="backward")


##############################################################################/
#Final model####
##############################################################################/

#Here are the results for the final model with the nAChR.81 genotype, 
#the CY3 expression level and their interaction
mod1<-glm(LC50~nAChR.81*CY3_EXP,data=sumDat,
          family=stats::gaussian(link="log"))
summary(mod1)
#anova(mod1,test="Chisq")
Anova(mod1,type=c("III"))
#plot of residuals versus fitted
plot(mod1,1)
#back transformation of the estimated coefficient
exp(coef(mod1))
exp(coef(mod1)+1.96*summary(mod1)$coefficients[,2])
exp(coef(mod1)-1.96*summary(mod1)$coefficients[,2])


##############################################################################/
#END
##############################################################################/