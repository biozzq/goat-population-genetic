#### get bin map
python3 02.CountReadsRatio.py -v divergence.vcf.gz -o Hybrid.CountReadsRatio
awk '$4>10' Hybrid.CountReadsRatio > Hybrid.CountReadsRatio.10
python3 04.ConvertToGT.py -i Hybrid.CountReadsRatio.10 -o Hybrid.CountReadsRatio.10.gt
sed s/X1/30/g Hybrid.CountReadsRatio.10.gt |sed s/X2/31/g|sort -k1,1n -k2,2n >Hybrid.CountReadsRatio.10.gt.sort
python3 06.MergeNonrecombinationRegion.py -i Hybrid.CountReadsRatio.10.gt.sort -o Hybrid.CountReadsRatio.10.gt.sort.merge
python3 07.PlotHeatmap.py -i Hybrid.CountReadsRatio.10.gt.sort -o Hybrid.CountReadsRatio.10.gt.sort.pdf

#### load the qtl package
rm(list=ls())
library(qtl)
Path="work path"       ## Work path
setwd(Path)
########## estimate the genetic map
mapthis <- read.cross("csv", ".", "PrimeryInput.csv",sep=",",estimate.map=T,map.function="kosambi",crosstype="f2")
write.cross(mapthis, format="csv",filestem="PrimeryInput-GM",digits=T)
########## change the genotype
mydata <- read.table("PrimeryInput-GM.csv", sep = ',',header=T)
mydata <- data.frame(lapply(mydata, function(x) {gsub("AA", "A", x)}))
mydata <- data.frame(lapply(mydata, function(x) {gsub("AB", "H", x)}))
mydata <- data.frame(lapply(mydata, function(x) {gsub("BB", "B", x)}))
mydata <- sapply(mydata, as.character)
mydata[is.na(mydata)] <- " "
write.table(mydata, "PrimeryInput-GM-edit.csv", sep = ',', row.names = F, col.names = T, quote = F)
########## reload the file with genetic map
mapthis <- read.cross("csv", ".", "PrimeryInput-GM-edit.csv")
summary(mapthis)
########## caculate the LOD value and plot
mapthis<-calc.genoprob(mapthis,step=5,error.prob = 1.0e-4,map.function = "kosambi",stepwidth = "fixed")
out<-scanone(mapthis,pheno.col=1,method="em")
operm<-scanone(mapthis,method = "em",n.perm = 500)
add.threshold(out,alpha = 0.05,perms = operm)
write.table(out, "Hybird.lod", sep = '\t', row.names = T, col.names = T, quote = F)
summary(out,perm=operm,alpha=0.05,format = "tabByChr",ci.function = "bayesint")
pdf("Hybird.lod.pdf")
plot(out)
dev.off()
