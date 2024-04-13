##############################################################################/
##############################################################################/
#Figure 4: proportion of metabolic resistance as a function of CYP6CY3
##############################################################################/
##############################################################################/

#first we need to run the regression analysis
source("MyMetab_reg_NTSR.R")


##############################################################################/
#Figure 4: Plotting the asymptotic regression model####
##############################################################################/

colovec=c("red3","orange3","green3")

pdf(file="output/Figure_4_metaRR.pdf",width=7,height=6)
op<-par(mar=c(5.1,5.1,1.1,1.1))
plot(modplot,type="confidence",log="",col="skyblue3",
     lwd=4,lty=2,ann=FALSE,axes=FALSE,bty="n",
     xlab="",
     ylab="",
     ylim=c(0,1),xlim=c(0,24))
points(sumDat$propMeta~sumDat$CY3_EXP,
       bg=colovec[as.numeric(sumDat$nAChR.81)],pch=21)
plot(modplot,type="obs",pch=19,col="green3",cex=1.5,add=TRUE)
axis(1,lwd=3,font=2)
axis(2,lwd=3,las=1,font=2)
box(bty="o",lwd=3)
title(xlab="CYP6CY3 expression level",
      ylab="Proportion of metabolic resistance",
      cex.lab=1.4,font=2,line=3.5)
legend(18,0.5,pch=c(21,21,19),pt.cex=c(1.5,1.5,2),
       pt.bg=colovec,col=c("black","black","green3"),
       legend=c(expression("81"^"TT"),
                expression("81"^"RT"),
                expression("81"^"RR")),
       y.intersp=1.4,bty="n")
par(op)
dev.off()


##############################################################################/
#END
##############################################################################/