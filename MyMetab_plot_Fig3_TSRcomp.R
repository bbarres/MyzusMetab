##############################################################################/
##############################################################################/
#Figure 3: comparison of TSR groups
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
source("MyMetab_load.R")


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