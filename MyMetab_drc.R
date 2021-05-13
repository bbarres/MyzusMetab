##############################################################################/
##############################################################################/
#Dose response curve analyses
##############################################################################/
##############################################################################/


#loading the packages necessary for the analysis
source("MyMetab_load.R")


##############################################################################/
#comparing the LD50 with or without PBO####
##############################################################################/

REZ<-data.frame("ech_id"=as.character(),"LD50"=as.character())

for (i in 1: dim(table(dataMyMeta$ech_id))[1]) {
  datatemp<-dataMyMeta[dataMyMeta$ech_id==names(table(dataMyMeta$ech_id))[i],]
  datatempOut<-datatemp[datatemp$synerg_id!="PBO",]
  temp.m1<-drm(nb_mtot/(nb_mtot+nb_vi)~dose,
               data=datatemp,curveid=synerg_id,
               fct=LN.3u())
  EDcomp(temp.m1,c(0.5,0.5),type="absolute")
  temp<-ED(temp.m1,0.5,type="absolute")
  plot(temp.m1,ylim=c(0,1),col=c("black","red"),type="all",
       lty=1,main=names(table(dataMyMeta$ech_id))[i])
  plot(temp.m1,ylim=c(0,1),col=c("black","red"),pch=c(16,17),
       lty=1,main=names(table(dataMyMeta$ech_id))[i],add=TRUE)
  tempx<-data.frame("ech_id"=names(table(dataMyMeta$ech_id))[i],
                    "ED50"=temp[1])
  REZ<-rbind(REZ,tempx)
}











##############################################################################/
#END
##############################################################################/