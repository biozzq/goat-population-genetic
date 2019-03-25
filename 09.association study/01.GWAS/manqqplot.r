source("~/script/GWAS/MANHATTAN_QQ.r")
setwd("./")

args<-commandArgs(TRUE)
IN<-args[1]
name<-args[2]

mydata=read.table(IN, header=TRUE)
GAPIT.Manhattan(mydata,name.of.trait=paste("",name,sep =""))
mydata=mydata[order(as.numeric(mydata[,3]),decreasing = FALSE),]
myoutdata=mydata[1:2001,]
myoutdata[,3] <-  -log10(myoutdata[,3])
write.table(myoutdata,file=paste("gwas.",name,".assoc.t2000.txt",sep =""), col.names=T, row.names=F, quote=F,sep=" ")

