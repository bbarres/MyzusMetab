##############################################################################/
##############################################################################/
#Dose response curve analyses
##############################################################################/
##############################################################################/


#loading the packages necessary for the analysis
source("MyMetab_load.R")


##############################################################################/
#Checking the quality of the control's bioassay for each date####
##############################################################################/

#first we limit the dataset to the control clone
ContCheData<-dataMyMeta[dataMyMeta$ech_id=="11-0037-0001",]
REZcont<-data.frame("ech_id"=as.character(),
                    "LD50"=as.character(),
                    "SE"=as.character(),
                    "LD50.PBO"=as.character(),
                    "SE.PBO"=as.character())

for (i in 1: dim(table(ContCheData$dat_test))[1]) {
  
  datatemp<-ContCheData[ContCheData$dat_test==
                          names(table(ContCheData$dat_test))[i],]
  temp.m1<-drm(nb_mtot/(nb_mtot+nb_vi)~dose,
               data=datatemp,curveid=synerg_id,
               fct=LN.3u())
  EDcomp(temp.m1,c(0.5,0.5),type="absolute")
  temp<-ED(temp.m1,0.5,type="absolute")
  plot(temp.m1,ylim=c(0,1),col=c("black","red"),type="all",
       lty=1,main=names(table(ContCheData$dat_test))[i],
       legendPos=c(50,1),ylab="Death rate")
  plot(temp.m1,ylim=c(0,1),col=c("black","red"),pch=c(16,17),
       lty=1,main=names(table(ContCheData$dat_test)[i]),
       legend=FALSE,add=TRUE)
  tempx<-data.frame("ech_id"=names(table(ContCheData$dat_test))[i],
                    "LD50"=temp[1],"SE"=temp[3],
                    "LD50.PBO"=temp[2],"SE.PBO"=temp[4])
  REZcont<-rbind(REZcont,tempx)

  }

plot(log(REZcont$LD50),ylab="Death rate")
plot(log(REZcont$LD50.PBO),ylab="Death rate")

dataMyCor<-dataMyMeta
#to remove : 20170503, 20170312 ?, 20180320
dataMyCor<-dataMyMeta[!dataMyMeta$dat_test %in% 
                        c("20170503","20170312","20180320") & 
                        dataMyMeta$ech_id!="11-0037-0001",]


##############################################################################/
#comparing the LD50 with or without PBO####
##############################################################################/

#removing experiment with more than 20% of mortality in the control dose
dataMyCor$mortrate<-dataMyCor$nb_mtot/(dataMyCor$nb_mtot+dataMyCor$nb_vi)
dataMyCor$idd<-paste(dataMyCor$ech_id,dataMyCor$dat_test,
                     dataMyCor$synerg_id,dataMyCor$rep_test,
                     sep=".")
repetorem<-dataMyCor[dataMyCor$mortrate>0.2 & dataMyCor$dose==0,]
dataMyCor<-dataMyCor[!dataMyCor$idd %in% repetorem$idd,]


REZ<-data.frame("ech_id"=as.character(),
                "LD50"=as.character(),
                "SE"=as.character(),
                "LD50.PBO"=as.character(),
                "SE.PBO"=as.character(),
                "EDcomp.est"=as.character(),
                "EDcomp.pval"=as.character())

for (i in 1: dim(table(dataMyCor$ech_id))[1]) {
  
  datatemp<-dataMyCor[dataMyCor$ech_id==names(table(dataMyCor$ech_id))[i],]
  datatempOut<-datatemp[datatemp$synerg_id!="PBO",]
  temp.m1<-drm(nb_mtot/(nb_mtot+nb_vi)~dose,
               data=datatemp,curveid=synerg_id,
               fct=LN.3u())
  temp0<-EDcomp(temp.m1,c(0.5,0.5),type="absolute")
  temp<-ED(temp.m1,0.5,type="absolute")
  plot(temp.m1,ylim=c(0,1),col=c("black","red"),type="all",
       lty=1,main=names(table(dataMyCor$ech_id))[i],
       legendPos=c(50,1),ylab="Death rate")
  plot(temp.m1,ylim=c(0,1),col=c("black","red"),pch=c(16,17),
       lty=1,main=names(table(dataMyCor$ech_id))[i],
       legend=FALSE,add=TRUE)
  tempx<-data.frame("ech_id"=names(table(dataMyCor$ech_id))[i],
                    "LD50"=temp[1],"SE"=temp[3],
                    "LD50.PBO"=temp[2],"SE.PBO"=temp[4],
                    "EDcomp.est"=temp0[1],
                    "EDcomp.pval"=temp0[4])
  REZ<-rbind(REZ,tempx)
  
}

#export the result table
write.table(REZ, file="output/results_bioassay.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


##############################################################################/
#END
##############################################################################/