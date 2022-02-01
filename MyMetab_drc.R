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

#first we limit the data set to the control clone
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

#to remove : 20170503, 20170312 ?, 20180320
dataMyCor<-dataMyMeta[!dataMyMeta$dat_test %in% 
                        c("20170503","20170312","20180320"),]


##############################################################################/
#comparing the LD50 with or without PBO####
##############################################################################/

# #if you want the clean data set by removing problematic date
# dataMyCor<-dataMyMeta[!dataMyMeta$dat_test %in% 
#                         c("20170503","20170312","20180320"),]
dataMyCor<-dataMyMeta
#removing experiment with more than 20% of mortality in the control dose
dataMyCor$mortrate<-dataMyCor$nb_mtot/(dataMyCor$nb_mtot+dataMyCor$nb_vi)
dataMyCor$effectif<-(dataMyCor$nb_mtot+dataMyCor$nb_vi)
dataMyCor$idd<-paste(dataMyCor$ech_id,dataMyCor$dat_test,
                     dataMyCor$synerg_id,dataMyCor$rep_test,
                     sep=".")
repetorem<-dataMyCor[dataMyCor$mortrate>0.2 & dataMyCor$dose==0,]
dataMyCor<-dataMyCor[!dataMyCor$idd %in% repetorem$idd,]


REZ<-data.frame("ech_id"=as.character(),
                "LD50"=as.character(),
                "SE"=as.character(),
                "mean/dose"=as.character(),
                "LD50.PBO"=as.character(),
                "SE.PBO"=as.character(),
                "mean/dose.PBO"=as.character(),
                "EDcomp.est"=as.character(),
                "EDcomp.pval"=as.character())

pdf(file="output/figure_by_clone.pdf",height=8,width=6)
nf<-layout(matrix(c(1,1,1,1,2,3),3,2,byrow=TRUE))
for (i in 1: dim(table(dataMyCor$ech_id))[1]) {
  
  datatemp<-dataMyCor[dataMyCor$ech_id==names(table(dataMyCor$ech_id))[i],]
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
  dattempWout<-datatemp[datatemp$synerg_id!="PBO",]
  temp1<-dattempWout[,c("ech_id","synerg_id","rep_test","test_echec",
                       "dose","effectif")]
  temp2<-pivot_wider(temp1,names_from="dose",values_from="effectif",
                     values_fn=mean)
  barplot(as.matrix(temp2[,c(5:dim(temp2)[2])]),las=2)
  dattempWith<-datatemp[datatemp$synerg_id=="PBO",]
  temp3<-dattempWith[,c("ech_id","synerg_id","rep_test","test_echec",
                       "dose","effectif")]
  temp4<-pivot_wider(temp3,names_from="dose",values_from="effectif",
                     values_fn=mean)
  barplot(as.matrix(temp4[,c(5:dim(temp4)[2])]),
          col=brewer.pal(9,"Reds")[c(7,5,3)],las=2)
  tempx<-data.frame("ech_id"=names(table(dataMyCor$ech_id))[i],
                    "LD50"=temp[1],"SE"=temp[3],
                    "mean/dose"=mean(colSums(
                      as.matrix(temp2[,c(5:dim(temp2)[2])])),
                      na.rm=TRUE),
                    "LD50.PBO"=temp[2],"SE.PBO"=temp[4],
                    "mean/dose.PBO"=mean(colSums(
                      as.matrix(temp4[,c(5:dim(temp4)[2])])),
                      na.rm=TRUE),
                    "EDcomp.est"=temp0[1],
                    "EDcomp.pval"=temp0[4])
  REZ<-rbind(REZ,tempx)
  
}

dev.off()

#export the result table
write.table(REZ, file="output/results_bioassay.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


##############################################################################/
#END
##############################################################################/