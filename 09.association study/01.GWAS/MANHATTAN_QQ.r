#Object: Manhattan & QQ Plot for genome-wide association study
#Output: Single figure with PDF format
#Authors: Zhiwu Zhang
#Modified by Zhengkui Zhou (zkzhou@126.com) 
#Last update: Feb 8, 2014




`GAPIT.Pruning` <-
function(values,DPP=5000){

if(length(values)<=DPP)return(c(1:length(values)))
  
#values= log.P.values
values=sqrt(values)  #This shift the weight a little bit to the low building.

#Handler of bias plot
rv=runif(length(values))
values=values+rv
values=values[order(values,decreasing = T)]

theMin=min(values)
theMax=max(values)
range=theMax-theMin
interval=range/DPP

ladder=round(values/interval)
ladder2=c(ladder[-1],0)
keep=ladder-ladder2
index=which(keep>0)


return(index)
}#end of GAPIT.Pruning 
#=============================================================================================


#Prune most non important SNPs off the plots, 
#Object: To get index of subset that evenly distribute
#`GAPIT.Pruning` <-  
#function(values,DPP){
#if(length(values)<=DPP)return(c(1:length(values)))
#values=sqrt(values)
#theMin=min(values)
#theMax=max(values)
#range=theMax-theMin
#interval=range/DPP
#ladder=round(values/interval)
#ladder2=c(ladder[-1],0)
#keep=ladder-ladder2
#index=which(keep>0)
#return(index)
#}

`GAPIT.Manhattan` <-
function(GI.MP = NULL, name.of.trait = "Trait",plot.type = "Genomewise", plot.type2 = "Chromosomewise", DPP=50000, cutOff=0.01, band=5, seqQTN=NULL){
print(paste("Start to plot figure for trait: ", name.of.trait ,sep = ""))
if(is.null(GI.MP)) return
P.values <- as.numeric(GI.MP[,3])
borrowSlot=4
GI.MP[,borrowSlot]=0 #Inicial as 0
if(!is.null(seqQTN))GI.MP[seqQTN,borrowSlot]=1

index=which(GI.MP[,borrowSlot]==1  & is.na(GI.MP[,3]))
GI.MP[index,3]=1
GI.MP=matrix(as.numeric(as.matrix(GI.MP) ) ,nrow(GI.MP),ncol(GI.MP))

#Remove all SNPs that do not have a choromosome, bp position and p value(NA)
GI.MP <- GI.MP[!is.na(GI.MP[,1]),]
GI.MP <- GI.MP[!is.na(GI.MP[,2]),]
GI.MP <- GI.MP[!is.na(GI.MP[,3]),]

#Remove all SNPs that have P values between 0 and 1 (not na etc)
GI.MP <- GI.MP[GI.MP[,3]>0,]
GI.MP <- GI.MP[GI.MP[,3]<=1,]

GI.MP <- GI.MP[GI.MP[,1]!=0,]

numMarker=nrow(GI.MP)
bonferroniCutOff=-log10(cutOff/numMarker)

#Replace P the -log10 of the P-values
GI.MP[,3] <-  -log10(GI.MP[,3])

y.lim <- as.integer(ceiling(max(GI.MP[,3])))
#y.lim = y.lim+1
print("The max -logP vlaue is")
print(y.lim)

chm.to.analyze <- unique(GI.MP[,1])
chm.to.analyze=chm.to.analyze[order(chm.to.analyze)]
numCHR= length(chm.to.analyze)

if(plot.type == "Genomewise")
{
nchr=max(chm.to.analyze)
ncycle=ceiling(nchr/band)
ncolor=band*ncycle
palette(rainbow(ncolor+1))
cycle1=seq(1,nchr,by= ncycle)
thecolor=cycle1

for(i in 2:ncycle){thecolor=c(thecolor,cycle1+(i-1))}
GI.MP <- GI.MP[order(GI.MP[,2]),]
GI.MP <- GI.MP[order(GI.MP[,1]),]
color.vector <- rep(c("orangered","navyblue"),numCHR)
ticks=NULL
lastbase=0

#change base position to accumulatives (ticks)
for (i in chm.to.analyze)
{
  index=(GI.MP[,1]==i)
  ticks <- c(ticks, lastbase+mean(GI.MP[index,2]))
  GI.MP[index,2]=GI.MP[index,2]+lastbase
  lastbase=max(GI.MP[index,2])
}
x0 <- as.numeric(GI.MP[,2]) ## BP
y0 <- as.numeric(GI.MP[,3]) ## P
z0 <- as.numeric(GI.MP[,1]) ## CHR
position=order(y0,decreasing = TRUE)
#print(length(position))
index0=GAPIT.Pruning(y0[position],DPP=DPP)
index=position[index0]
x=x0[index]
y=y0[index]
z=z0[index]
output=cbind(z,x,y)
#print(output)
#write.table(output, paste("gwasplot.","printy",".txt",sep=""), sep = '\t', row.names = F, col.names = T, quote = F)
#Extract QTN
QTN=GI.MP[which(GI.MP[,borrowSlot]==1),]

#Draw circles with same size and different thikness
size=1
ratio=5
base=1
themax=max(y)
themin=min(y)
wd=((y-themin+base)/(themax-themin+base))*size*ratio
s=size-wd/ratio/2
y.lim <- as.integer(ceiling(max(GI.MP[,3])))
#y.lim = 10
  pdf(paste("GWAS.", name.of.trait,".Manhattan-QQ Plot.pdf" ,sep = ""), width = 18,height=4.5)
  layout(matrix(c(1,2), 1, 2, byrow = TRUE),widths=c(3,1), heights=c(1,1))  ##  
  par(mar = (c(3,4,2,2)+ 0.5), mgp=c(1.6,1,0))   ##  
  par(bty="l", lwd=1.5)  ## bty=l  the plot is coordinate instead of box
  mycols=rep(c("orangered","navyblue"),max(z))  ## doulb dolor loop by chromosome 
plot(y~x,ylab=expression(-log[10](italic(p))) , ylim=c(0,y.lim), xaxs="i", yaxs="i" ,
     cex.axis=0.5, cex.lab=1.0 ,col=mycols[z], axes=FALSE, type = "p",
     pch=20,lwd=wd,cex=0.5, xlab="Chromosome")
if(!is.null(dim(QTN)))abline(v=QTN[,2], lty = 2, lwd=1.5, col = "grey")
abline(h=bonferroniCutOff,col="dimgray",lty=2, lwd=1)
title(xlab="Chromosome")
 axis(1, at=ticks,tck=-0.01, cex.axis=1,labels=chm.to.analyze,tick=T, lwd=1.5, padj=-1)    ## lwd: line width tick=T, 
 axis(2, tck=-0.01, cex.axis=1,lwd=1.5, padj=1)     ## tck=-0.01 let the tck shor            
 box()
 mtext(paste("GWAS on ",name.of.trait, sep=""), cex=1.5, font.main =1,  side=3, outer=TRUE,line=-1.5) 
print("Manhattan-Plot.Genomewise finished!")

############# QQ-plot#################
DPP=50000
if(length(P.values[P.values>0])<1) return(NULL)
DPP=round(DPP/4) #Reduce to 1/4 for QQ plot
P.values <- P.values[order(P.values)]
p_value_quantiles <- (1:length(P.values))/(length(P.values)+1)
    log.P.values <- -log10(P.values)
    log.Quantiles <- -log10(p_value_quantiles)
    index=GAPIT.Pruning(log.P.values,DPP=DPP)
    log.P.values=log.P.values[index ]
    log.Quantiles=log.Quantiles[index]
    qqplot(log.Quantiles, log.P.values, xlim = c(0,max(log.Quantiles)), ylim = c(0,max(log.P.values)),pch=20, mgp=c(1.6,0.5,0),
           cex.axis=1, cex.lab=1, lty = 1, lwd = 1.2, col = "Blue" ,xlab =expression(Expected~~-log[10](italic(p))), tck=-0.015,
           ylab = expression(Observed~~-log[10](italic(p)))) 
    abline(a = 0, b = 1, col = "red")
dev.off()
print("QQ-Plot finished!")
print(paste("Manhattan & QQ Plot for trait: ", name.of.trait," accomplished successfully!" ,sep = ""))
}
}
