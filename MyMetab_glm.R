##############################################################################/
##############################################################################/
#Modelling the LD50
##############################################################################/
##############################################################################/


#loading the packages necessary for the analysis
source("MyMetab_load.R")


##############################################################################/
#Regression analyses####
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
#export to .pdf 9 x 7 inches

#because the correlation between copy number of the different genes
#and level of expression of the different genes are correlated, we 
#pick only one of these variables. Since CY23 is a less good candidate
#we don't include the related variables either

mod1<-glm(LC50~genetic.group+nAChR.81*CY3_CN,data=sumDat,
          family=stats::gaussian(link="log"))
summary(mod1)
#the genetic.group is not important, so we remove this 
#variable from the model

mod1<-glm(LC50~nAChR.81*CY3_EXP,data=sumDat,
          family=stats::gaussian(link="log"))
summary(mod1) #best model AIC 320
mod2<-glm(LC50~nAChR.81*CY4_EXP,data=sumDat,
          family=stats::gaussian(link="log"))
summary(mod2) #AIC 390
mod3<-glm(LC50~nAChR.81+CY3_EXP,data=sumDat,
          family=stats::gaussian(link="log"))
summary(mod3) #AIC 351
mod4<-glm(LC50~nAChR.81+CY4_EXP,data=sumDat,
          family=stats::gaussian(link="log"))
summary(mod4) #AIC 389

#modelling the EC50 with PBO
modPBO.1<-glm(LC50.PBO~nAChR.81*CY3_EXP,data=sumDat,
          family=stats::gaussian(link="log"))
summary(modPBO.1) #best model AIC 282
modPBO.2<-glm(LC50.PBO~nAChR.81*CY4_EXP,data=sumDat,
          family=stats::gaussian(link="log"))
summary(modPBO.2) #AIC 324
modPBO.3<-glm(LC50.PBO~nAChR.81+CY3_EXP,data=sumDat,
          family=stats::gaussian(link="log"))
summary(modPBO.3) #AIC 292
modPBO.4<-glm(LC50.PBO~nAChR.81+CY4_EXP,data=sumDat,
          family=stats::gaussian(link="log"))
summary(modPBO.4) #AIC 322

#modelling the difference/ratio with or without PBO
modDif1<-glm(LC50/LC50.PBO~nAChR.81*CY3_EXP,data=sumDat,
             family=stats::gaussian(link="identity"))
summary(modDif1) #best model AIC 76
modDif2<-glm(LC50/LC50.PBO~nAChR.81*CY4_EXP,data=sumDat,
             family=stats::gaussian(link="identity"))
summary(modDif2) #AIC 77
modDif3<-glm(LC50/LC50.PBO~nAChR.81+CY3_EXP,data=sumDat,
             family=stats::gaussian(link="identity"))
summary(modDif3) #AIC 80 
modDif4<-glm(LC50/LC50.PBO~nAChR.81+CY4_EXP,data=sumDat,
             family=stats::gaussian(link="identity"))
summary(modDif4) #AIC 81

#to sum up the best model
summary(mod1)
summary(modPBO.1)
summary(modDif1)


##############################################################################/
#END
##############################################################################/