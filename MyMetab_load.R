##############################################################################/
##############################################################################/
#Data and package loading for the analyses and figures plotting
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
library(drc)
library(medrc)
library(plotrix)
library(gdata)
library(tidyr)
library(RColorBrewer)


##############################################################################/
#loading the bioassay data set####
##############################################################################/

#load the global data set
dataMyMeta<-read.table(file="data/donnees_Myzus_P450_20201224_2.txt",
                       header=T,sep=";")

#because some concentration were only used for adapting the pesticide dose
#scale. It also include repetition that were flawed because of insufficient
#number of individual for the entire repetition or because there was a 
#problem during the lab experiment
dataMyMeta<-dataMyMeta[dataMyMeta$test_echec!=1,]



# Pour l’analyse, nous avions décidé des critères suivants :
#   
#   -          “For each concentration that was tested, three replicates or more, 
# involving at least 10 L1 per replicate, were tested. Tests were repeated until
# a minimum of 45 aphids per dose were assayed for each modality.”
# 
# -          Mortalité dans le témoin < 20%
# 
# -          Pour qu’un test sur un clone donné soit valide, le test sur la 
# référence 11-0037-0001 correspondant doit être valable.
# 
# -          En plus de ces critères stricts, lors des premières analyses R, 
# nous supprimions parfois les données des doses extrêmes lorsque le plateau 
# était trop grand, pour que le modèle gère mieux ces résultats. Nous faisions 
# ça manuellement.



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