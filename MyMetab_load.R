##############################################################################/
##############################################################################/
#Data and package loading for the analyses and figures plotting
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
library(car)
library(drc)
library(emmeans)
library(gdata)
library(jtools)
library(lattice)
library(medrc)
library(plotrix)
library(RColorBrewer)
library(tidyr)
library(vioplot)


##############################################################################/
#loading the bioassay data set####
##############################################################################/

#load the global data set
dataMyMeta<-read.table(file="data/bioassayRawData.txt",
                       header=TRUE,sep=";")

#because some concentration were only used for adapting the pesticide dose
#scale. It also include repetition that were flawed because of insufficient
#number of individual for the entire repetition or because there was a 
#problem during the lab experiment
dataMyMeta<-dataMyMeta[dataMyMeta$test_echec!=1,]

# For the cleaning of the data, we used these criteria :
# 
# -“For each concentration that was tested, three replicates or more,
# involving at least 10 L1 per replicate, were tested. Tests were repeated 
# until a minimum of 45 aphids per dose were assayed for each modality.”
# 
# -untreated control mortality < 20%
# 
# -to validate a test on a clone, the reference clone 11-0037-0001 tested 
# at the same date should also be valid

#load data for the regression model
sumDat<-read.table(file="data/summaData.txt",header=TRUE,sep="\t",
                   stringsAsFactors=TRUE)
sumDat$nAChR.81<-factor(sumDat$nAChR.81,levels=c("TT","RT","RR"))
levels(sumDat$nAChR.81)<-c("[TT]","[RT]","[RR]")
#adding columns for comparisons with and without PBO
sumDat$propMeta<-(1-(sumDat$LC50.PBO/sumDat$LC50))
sumDat$diffMeta<-(sumDat$LC50-sumDat$LC50.PBO)


##############################################################################/
#Writing info session for reproducibility####
##############################################################################/

sink("session_info.txt")
print(sessioninfo::session_info())
sink()
#inspired by an R gist of François Briatte: 
#https://gist.github.com/briatte/14e47fb0cfb8801f25c889edea3fcd9b


##############################################################################/
#END
##############################################################################/