##############################################################################/
##############################################################################/
#Comparison between LD50 with or without PBO 
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
source("MyMetab_load.R")


##############################################################################/
#Correlation between LD50 with or without PBO####
##############################################################################/

#one way anova
simp.mod<-lm(LC50.PBO~nAChR.81,data=sumDat)
summary(simp.mod)
#post hoc tests
TukeyHSD(aov(simp.mod))
plot(simp.mod,1) #seems there is different variance between group
#testing heteroscedasticity
leveneTest(LC50.PBO~nAChR.81,data=sumDat) #there is heteroscedasticity

#log transformation to correct for heteroscedasticity
log.mod<-lm(log(LC50.PBO)~nAChR.81,data=sumDat)
summary(log.mod)
#post hoc tests
TukeyHSD(aov(log.mod))
plot(log.mod,1)
#Means and SE obtained by back-transformation via the delta method
BackTrans<-emmeans(log.mod,~nAChR.81,type="response")

#anova with Welch correction
oneway.test(sumDat$LC50.PBO~sumDat$nAChR.81,var.equal=FALSE)
pairwise.t.test(sumDat$LC50.PBO,sumDat$nAChR.81,
                pool.sd=FALSE,p.adjust.method="BH")

#another more conservative methods: non parametric test
kruskal.test((sumDat$LC50.PBO)~sumDat$nAChR.81)
#pairwise test
pairwise.wilcox.test(sumDat$LC50.PBO,sumDat$nAChR.81,p.adjust.method="BH")


##############################################################################/
#LC50 dumbellplot with or without PBO####
##############################################################################/

#this function was retrieved from:
#https://r-coder.com/dot-plot-r/?utm_content=cmp-true

# v1: numeric variable
# v2: numeric variable
# group: vector (numeric or character) or a factor containing groups
# labels: labels for the dot chart
# segments: whether to add segments (TRUE) or not (FALSE)
# text: whether to add text (TRUE) or not (FALSE)
# pch: symbol
# col1: color of the variable v1. If you want to
# add group colors add them here
# col1: color of the variable v2
# pt.cex: size of the points
# segcol: color of the segment
# lwd: width of the segment
# ... : additional arguments to be passed to dotchart function

dumbbell <- function(v1, v2, group = rep(1, length(v1)), labels = NULL,
                     segments = FALSE, text = FALSE, pch = 19,
                     colv1 = 1, colv2 = 1, pt.cex = 1, segcol = 1,
                     lwd = 1,logdum="", ...) {
  
  o <- sort.list(as.numeric(group), decreasing = TRUE)
  group <- group[o]
  offset <- cumsum(c(0, diff(as.numeric(group)) != 0))
  y <- 1L:length(v1) + 2 * offset
  
  dotchart(v1, labels = labels, color = colv1, xlim = range(v1, v2) + c(-2, 2),
           groups = group, pch = pch, pt.cex = pt.cex,log=logdum)
  
  if(segments == TRUE) {
    for(i in 1:length(v1)) {
      segments(min(v2[i], v1[i]), y[i],
               max(v2[i], v1[i]), y[i],
               lwd = lwd, col = segcol) 
    }
  }
  
  for(i in 1:length(v1)){
    points(v2[i], y[i], pch = pch, cex = pt.cex, col = colv2)
    points(v1[i], y[i], pch = pch, cex = pt.cex, col = colv1)
  }
  
  if(text == TRUE) {
    for(i in 1:length(v1)) {
      text(min(v2[i ], v1[i]) - 1.5, y[i],
           labels = min(v2[i], v1[i]))
      text(max(v2[i], v1[i]) + 1.5, y[i],
           labels = max(v2[i], v1[i])) 
    }
  }
}



x <- sumDat[order(sumDat$LC50), ] 

dumbbell(v1 = x$LC50.PBO, v2 = x$LC50, group = x$nAChR.81,
         text = TRUE, segcol = "gray", lwd = 3, labels = x$clone.ID,
         segments = TRUE, pch = 19, pt.cex = 1.5, colv1 = 1, colv2 = "blue",
         logdum="x")


#pdf(file="output/Figure_X_TSRcomp.pdf",width=6,height=5.5)
#op<-par(mar=c(5.1,5.1,1.1,1.1))
col=c("green3","orange3","red3")
c("green4","orange4","red4"))

dotchart(sort(sumDat$LC50),
         groups=sumDat$nAChR.81[order(sumDat$LC50)],
         labels=sumDat$clone.ID[order(sumDat$LC50)],log="x")
dotchart(sumDat$LC50.PBO[order(sumDat$LC50)],add=TRUE)

axis(1,at=c(1,2,3),labels=c("[RR]","[RT]","[TT]"),lwd=3,font=2)
axis(2,at=c(log(50),log(100),log(200),log(500),log(1000),log(2000)),
     labels=c("50","100","200","500","1000","2000"),lwd=3,
     las=1,font=2)
title(xlab="Target-site resistance genotype",
      ylab="LC50 with PBO",
      cex.lab=1.5,font=2)
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